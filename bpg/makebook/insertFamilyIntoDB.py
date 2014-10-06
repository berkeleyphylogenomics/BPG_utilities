#!/usr/bin/env python

import cPickle
import datetime
from optparse import OptionParser
import os
import string
import sys
import re

from Bio import SeqIO, AlignIO
from Bio.SeqUtils import CheckSum
from django.db import connection, transaction
from psycopg2 import IntegrityError

from bpg.common.parsers.hmm_parsers \
    import parse_results_of_hmmsearch_or_hmmscan
from bpg.common.utils.compute_alignment_conservation \
    import update_tree_node_alignment_conservation
from pfacts003.phylofacts.models import Sequence, AlignedSequence,\
    SequenceHeader, UniProt, ThirdPartySequenceHeader, ThirdParty,\
    Family, BuildAlignmentNotes, Tree, TreeNode, TreeNodeAlignment,\
    TreeNodeSignalPeptide, HMM, HMM_Consensus, TreeNodeConsensus,\
    TreeNodeTransmembrane, UniProtDatIndex, Pfam, SequenceHMM, PDB_Chain,\
    TreeNodePDB
from pfacts003.utils.id_patterns import uniprot_accession_pat1,\
                                        uniprot_accession_pat2,\
                                        swissprot_desc_pat, uniprot_taxon_pat

current_pfam_version = 24.0

id_re = re.compile('SEQ[0-9]*')
seqid_re = re.compile('(SEQ[0-9]*)')
seqhdr_re = re.compile('(SEQHDR[0-9]*)')
capturing_seqhdr_re = re.compile('(SEQHDR([0-9]*))')
branch_length_re = re.compile('-?[0-9]+\.?[0-9e\-]*')
likelihood_value_re = re.compile('[0-9]\.[0-9e\-]*')
uppercase_translation = string.maketrans(string.lowercase, string.uppercase)
dotdash = '.-'
trivial_translation = string.maketrans('', '')
dotlowercase = string.lowercase + '.'
# This module's versions of these regular expressions don't require them to
# match the entire string
uniprot_accession_re1 = re.compile(uniprot_accession_pat1)
uniprot_accession_re2 = re.compile(uniprot_accession_pat2)
uniprot_identifier_re3 = re.compile('%s_%s' % (swissprot_desc_pat,
                                              uniprot_taxon_pat))


class node:
  def __init__(self, seqid=None):
    self.parent = None
    self.children = set()
    self.leftId = None
    self.rightId = None
    self.level = None
    self.branch_length = -1.0
    self.seqid = seqid
    self.tree_node = None
    self.likelihood = -1.0
    if seqid:
      self.contained_seqids = set([seqid])
    else:
      self.contained_seqids = set()

  def addChild(self, child):
    self.children.add(child)
    child.parent = self
    self.contained_seqids = self.contained_seqids | child.contained_seqids

  def updateLeftId(self, leftId, level):
    self.leftId = leftId
    self.level = level
    rightId = leftId + 1
    for child in self.children:
      rightId = child.updateLeftId(rightId, level + 1)
    self.rightId = rightId
    return rightId + 1

  def readFromTreeString(self, tree_string, i):
    if i >= len(tree_string):
      return len(tree_string)
    if tree_string[i] == '(':
      while tree_string[i] != ')':
        child = node()
        i = child.readFromTreeString(tree_string, i+1)
        self.addChild(child)
      i += 1
      m = likelihood_value_re.match(tree_string[i:])
      if m:
        self.likelihood = float(m.group(0))
        i += len(m.group(0))
    else:
      m = id_re.match(tree_string[i:])
      seqid = m.group(0)
      i += len(seqid)
      self.seqid = seqid
      self.contained_seqids = set([self.seqid])
    if i < len(tree_string) and tree_string[i] == ':':
      i += 1
      branch_length = branch_length_re.match(tree_string[i:]).group(0)
      self.branch_length = float(branch_length)
      i += len(branch_length)
    return i

  def printSelf(self, indentLevel = 0):
    if self.seqid:
      sys.stdout.write("%s" % self.seqid)
    else:
      sys.stdout.write("(")
      children = list(self.children)
      for i in range(len(children)):
        children[i].printSelf(indentLevel + 1)
        if i < len(children) - 1:
          sys.stdout.write(",")
      sys.stdout.write(")")
    if self.likelihood >= 0.0:
      sys.stdout.write("%g" % self.likelihood)
    if self.branch_length >= 0.0:
      sys.stdout.write(":%g" % self.branch_length)

  def createTreeNodeObjects(self, tree, sequence_header_of_seqid):
    parent_left_id = 0
    if self.parent:
      parent_left_id = self.parent.leftId
    sequence_header = None
    if self.seqid:
      sequence_header = sequence_header_of_seqid[self.seqid]
    self.tree_node = TreeNode.objects.create(tree=tree,
                                        left_id=self.leftId,
                                        right_id=self.rightId,
                                        parent_left_id=parent_left_id,
                                        level=self.level,
                                        sequence_header=sequence_header,
                                        distance_to_parent=self.branch_length)
    if self.likelihood >= 0.0:
      self.tree_node.likelihood = self.likelihood
      self.tree_node.save()
    for child in self.children:
      child.createTreeNodeObjects(tree, sequence_header_of_seqid)

  def getContainedLeaves(self):
    if self.seqid:
      return set([self])
    else:
      contained_leaves = set()
      for child in self.children:
        contained_leaves = contained_leaves | child.getContainedLeaves()
      return contained_leaves

@transaction.commit_on_success
def _sequence_header(header, sequence, uniprot=None, taxon=None):
    """Internal helper function to create sequence header

    When sharding this script, some code logic needed to be isolated
    into smaller executable pieces. This is so that, when running
    several instances of the same script on the same database tables,
    we don't get race conditions updating the same records.

    This routine takes a header (string), a sequence The SequenceHeader
    object has four main attributes (description, sequence, uniprot and
    taxon). Description and sequence are required and are given to this
    routine.
    """
    is_dirty = False

    sequence_header, is_created =\
        SequenceHeader.objects.get_or_create(header=header,
                                             sequence=sequence)

    if uniprot is not None:
        sequence_header.uniprot = uniprot
        is_dirty = True

    if taxon is not None:
        sequence_header.taxon = taxon
        is_dirty = True

    # Don't make expensive database call unless it is necessary
    if is_dirty:
        sequence_header.save()

    return sequence_header


@transaction.commit_manually
def _sequence(obj, chars, seguid):
    """Internal helper function to create sequence or aligned_sequence

    When sharding this script, some code logic needed to be isolated
    into smaller executable pieces. This is so that, when running
    several instances of the same script on the same database tables,
    we don't get race conditions updating the same records.

    This routine takes the object (obj) that is either a Sequence
    Django Model object or an AlignedSequence Django Model Object,
    an amino acid character string (chars) and a unique hash on
    those characters (seguid). All three parameters are required for
    this function.
    
    Debugging Notes: More care is needed with the get_or_create method
    than is generally true in other places. get_or_create works on the 
    primary key and should eliminate any sharding race conditions.
    However, the seguid is not the primary key in this case (obviously
    for performance reasons). So, there is a snapshot in time (from the
    get_or_create until the sequence.save() statements) where both
    sequence.chars and sequence.seguid will be the empty string.
    Because there cannot be two identical sequence.seguids (i.e., both
    empty string), we must lock the table for a small duration of time
    until the sequence.seguid gets updated.

    This is wrapped in a transaction that should roll back upon error.
    However, if this routine has been altered, and there are violations
    for unique constraints on seguid, check to verify there is not still
    an empty string for any seguid. 

    DO NOT make a check constraint against the empty string on
    sequence.seguid since the get_or_create method needs to have an
    empty field for this time for the short amount of time until the
    object is saved.
    """

    sid = transaction.savepoint()
    retry = True

    # This loop implements a tight spin-lock type construct
    # that should only spin twice at most. In practice this
    # has been shown to be effective on 20,000 PFAM families
    # during insertion.
    while retry is True:
        sequence, is_created =\
            obj.objects.get_or_create(seguid__exact=seguid)

        # The following line is NOT redundant. It *is* necessary.
        sequence.seguid = seguid
        sequence.chars = chars

        try:
            sequence.save()
        except IntegrityError:
            transaction.savepoint_rollback(sid)
            print "Retrying save of '%s' for '%s'" % (obj._meta.db_table,
                                                      seguid)
            retry = True
        else: 
            transaction.savepoint_commit(sid)
            retry = False

    transaction.commit()
    return sequence


def create_sequence_header(identifier, description, sequence):
  """Either find or create a Django SequenceHeader object

  Given the identifier, description and sequence, find or create a
  Sequence Header Django model object for the sequence given.  The
  identifier and description are from the FASTA sequence defline.
  Identifier is a parsed portion of description.  Sequence is a Django
  Sequence object.

  The SequenceHeader object has four main attributes (description,
  sequence, uniprot and taxon). Description and sequence are required
  and are given to this routine. Both UniProt and taxon come from a
  Django UniProt object.

  The appropriate UniProt object is searched for exhaustively to give
  the two optional attributes (uniprot and taxon). If the appropriate
  UniProt object is not found, a SequenceHeader object is created
  without these optional attributes.

  Because this application is ran in parallel, it is very possible to
  execute many instances of this routine at the exact moment in time.
  Therefore, one cannot assume that a search moments before is still
  valid. To mimimize this:
     * SequenceHeader.objects.get_or_create is used instead of first 
       getting the object and then later updating the object.
     * This logic was reduced to a single creation at the bottom of
       the routine
     * For simplicity, all UniProt data is assumed to be static. This
       is generally updated rarely in off-peak times. Note that any
       update to UniProt at the same time as this script could
       theoretically cause this script to fail.
  """
  # Try to find the UniProt record first

  uniprot_accession = ''
  # Check for a UniProt accession in the identifier or description
  m = uniprot_accession_re1.search(identifier)
  if m:
    uniprot_accession = m.group(0)
  else:
    m = uniprot_accession_re2.search(identifier)
    if m:
      uniprot_accession = m.group(0)
    else:
      m = uniprot_accession_re1.search(description)
      if m:
        uniprot_accession = m.group(0)
      else:
        m = uniprot_accession_re2.search(description)
        if m:
          uniprot_accession = m.group(0)

  if len(uniprot_accession) > 0:
    try:
      uniprot = UniProt.objects.get(accession__exact = uniprot_accession)
      sequence_header = _sequence_header(description, sequence,
                                    uniprot=uniprot, taxon=uniprot.taxon)

      return sequence_header
    except UniProt.DoesNotExist:
      try:
        # This must be an old accession, look for it in the uniprot_dat_index
        # table
        uniprot_dat_index = UniProtDatIndex.objects.get(uniprot_accession__exact
                                                            = uniprot_accession)
        uniprot = uniprot_dat_index.uniprot
        sequence_header = _sequence_header(description, sequence,
                                           uniprot=uniprot, taxon=uniprot.taxon)
        return sequence_header
      except UniProtDatIndex.DoesNotExist:
        # Maybe this accession is obsolete.  That didn't work, keep trying
        pass

  # Check for a SwissProt identifier
  uniprot_identifier = ''
  m = uniprot_identifier_re3.search(identifier)
  if m:
    uniprot_identifier = m.group(0)
  else:
    m = uniprot_identifier_re3.search(description)
    if m:
      uniprot_identifier = m.group(0)

  if len(uniprot_identifier) > 0:
    try:
      uniprot = UniProt.objects.get(uniprot_identifier__exact 
                                                          = uniprot_identifier)
      sequence_header = _sequence_header(description, sequence,
                                         uniprot=uniprot, taxon=uniprot.taxon)
      return sequence_header
    except UniProt.DoesNotExist:
      # That didn't work, keep trying
      pass

  # Look for Uniprot records with the same sequence
  seguid = sequence.seguid
  uniprot_objects = UniProt.objects.filter(seguid__exact = seguid)
  if uniprot_objects:
    if len(uniprot_objects) == 1:
      uniprot = uniprot_objects[0]
      sequence_header = _sequence_header(description, sequence,
                                         uniprot=uniprot, taxon=uniprot.taxon)
      return sequence_header

    # For each of the uniprot objects, try to find its taxon in the
    # description 
    taxa = set()
    max_matching_taxon_length = 0
    uniprot_with_max_matching_taxon_length = None
    for uniprot in uniprot_objects:
      taxa.add(uniprot.taxon)
      if uniprot.taxon.scientific_name and \
          description.find(uniprot.taxon.scientific_name):
        if len(uniprot.taxon.scientific_name) > max_matching_taxon_length:
          max_matching_taxon_length = len(uniprot.taxon.scientific_name)
          uniprot_with_max_matching_taxon_length = uniprot
      if uniprot.taxon.common_name and \
          description.find(uniprot.taxon.common_name):
        if len(uniprot.taxon.common_name) > max_matching_taxon_length:
          max_matching_taxon_length = len(uniprot.taxon.common_name)
          uniprot_with_max_matching_taxon_length = uniprot
    if max_matching_taxon_length > 0:
      uniprot = uniprot_with_max_matching_taxon_length
      sequence_header = _sequence_header(description, sequence,
                                         uniprot=uniprot, taxon=uniprot.taxon)
      return sequence_header
    elif len(taxa) == 1:
      # It's rather odd that there are multiple Uniprot entries with the same
      # sequence and the same taxon, but anyway, we can just pick one
      sequence_header = _sequence_header(description, sequence,
                                         uniprot=uniprot, taxon=uniprot.taxon)
      return sequence_header
    else:
      # There are multiple Uniprot entries with the same sequence and different
      # taxa, none of which is specified in the description of this sequence,
      # so unfortunately we can't identify either its Uniprot record or taxon
      sequence_header = _sequence_header(description, sequence)
      return sequence_header

  # Look for third-party records with the same sequence
  thirdparty_sequence_header_objects = ThirdPartySequenceHeader.objects.filter(
      sequence_header__sequence__exact = sequence)

  # See if all these third-party records come from a single taxon
  taxa = set()
  for thirdparty_sequence_header in thirdparty_sequence_header_objects:
    if thirdparty_sequence_header.sequence_header.taxon:
      taxa.add(thirdparty_sequence_header.sequence_header.taxon)

  if len(taxa) == 1:
    taxon = list(taxa)[0]
    sequence_header = _sequence_header(description, sequence,
                                       taxon=uniprot.taxon)

    return sequence_header

  # We've tried everything, and can't determine either the Uniprot record or
  # the taxon
  sequence_header = _sequence_header(description, sequence)
  return sequence_header


def insertPhobiusPredictionsIntoDB(canonical_root_node, phobius_filename):
  f = open(phobius_filename)
  lines = f.readlines()
  f.close()
  for line in lines:
    if line[0:2] == 'FT':
      fields = line.split()
      if len(fields) >= 4:
        if fields[1] == 'SIGNAL':
          start = int(fields[2])
          end = int(fields[3])
          TreeNodeSignalPeptide.objects.create(seq_start = start, 
                                              seq_end = end,
                                              method = 'phobius',
                                              tree_node = canonical_root_node)
        elif fields[1] == 'TRANSMEM':
          start = int(fields[2])
          end = int(fields[3])
          TreeNodeTransmembrane.objects.create(seq_start = start, 
                                              seq_end = end,
                                              method = 'phobius',
                                              tree_node = canonical_root_node)

def replace_seqids_by_seq_header_ids(in_tree_filename, out_tree_filename,
                                      sequence_header_of_seqid):
  f = open(in_tree_filename)
  tree_str = f.read()
  f.close()
  tree_tokens = seqid_re.split(tree_str)
  for i in range(len(tree_tokens)):
    if i % 2 == 1:
      tree_tokens[i] = "SEQHDR%d" % sequence_header_of_seqid[tree_tokens[i]].id
  new_tree_str = ''.join(tree_tokens)
  outf = open(out_tree_filename, "w")
  outf.write(new_tree_str)
  outf.close()

def insertPFAMPredictionsIntoDB(consensus_sequence, basename):
  # This file is one of the results of running hmmscan during buildFamily.
  hmmscan_filename = basename + "_vs_Pfam-A.hmmscan.out"
  hmmscan_results = parse_results_of_hmmsearch_or_hmmscan.parse(
                      hmmscan_filename, 0.001, 1, 1)
  # There should only be one query here - the family consensus sequence
  for query in hmmscan_results.hit_result_of_name_of_query:
    for pfam_name in hmmscan_results.hit_result_of_name_of_query[query]:
      hit_result = hmmscan_results.hit_result_of_name_of_query[query][pfam_name]
      pfam_objects = Pfam.objects.filter(name__exact = pfam_name,
                            overall_pfam_version__exact = current_pfam_version)
      if pfam_objects:
        hmm_objects = HMM.objects.filter(pfam__exact = pfam_objects[0])
        if hmm_objects:
          hmm_object = hmm_objects[0]
        else:
          print "No HMM entry for Pfam entry %s.  Reload Pfam" % pfam_name
          continue
        for match_number in hit_result.matches:
          match_result = hit_result.matches[match_number]
          SequenceHMM.objects.create(hmm = hmm_object,
                        sequence = consensus_sequence,
                        bit_score = match_result.bit_score,
                        e_value = match_result.i_evalue,
                        sequence_type = 'consensus sequence (different family)',
                        hmm_start = match_result.hmm_from,
                        hmm_end = match_result.hmm_to,
                        sequence_start = match_result.seq_from,
                        sequence_end = match_result.seq_to,
                        match_type = match_result.match_type,
                        n_aligned_chars = match_result.num_aligned_chars)
      else:
        print "Unrecognized Pfam name %s for Pfam version %0.1f." \
              % (pfam_name, current_pfam_version)
        print "Check that this is the current Pfam version",
        print "and whether the pfam table is up-to-date."

def insertPDBPredictionsIntoDB(hmm, tree_node, basename):
    hmmsearch_filename = basename + "_vs_PDB.hmmsearch.out"
    hmmsearch_results = parse_results_of_hmmsearch_or_hmmscan.parse(
                      hmmsearch_filename, 0.001, 1, 1)
    # There should only be one query here - the family HMM
    for query in hmmsearch_results.hit_result_of_name_of_query:
      for pdb_chain_id in hmmsearch_results.hit_result_of_name_of_query[query]:
        pdb_id, chain_id = pdb_chain_id.split('_')
        pdb_chain_objects = PDB_Chain.objects.filter(pdb__id__exact = pdb_id,
                                                    chain_id__exact = chain_id)
        if pdb_chain_objects:
          pdb_chain = pdb_chain_objects[0]
        else:
          print "Unrecognized PDB chain %s in hmmsearch results." \
            % pdb_chain_id,
          print "The PDB_Chain table may be out of date."
          continue
        hit_result \
            = hmmsearch_results.hit_result_of_name_of_query[query][pdb_chain_id]
        for match_number in hit_result.matches:
          match_result = hit_result.matches[match_number]
          aligned_seguid = CheckSum.seguid(match_result.aligned_hit)
          aligned_sequence_objects = AlignedSequence.objects.filter(
                                                seguid__exact = aligned_seguid)
          if aligned_sequence_objects:
            aligned_sequence = aligned_sequence_objects[0]
          else:
            # Because multiple versions of this fun run simultaneously, it is
            # possible this was just created moments ago. So, get_or_create
            # just in case.
            aligned_sequence, is_created  \
                = AlignedSequence.objects.get_or_create(
                                              chars = match_result.aligned_hit,
                                              seguid = aligned_seguid)
          sequence_hmm = SequenceHMM.objects.create(hmm = hmm,
                              sequence = pdb_chain.full_sequence,
                              aligned_sequence = aligned_sequence,
                              bit_score = match_result.bit_score,
                              e_value = match_result.i_evalue,
                              sequence_type = 'query',
                              hmm_start = match_result.hmm_from,
                              hmm_end = match_result.hmm_to,
                              sequence_start = match_result.seq_from,
                              sequence_end = match_result.seq_to,
                              match_type = match_result.match_type,
                              n_aligned_chars = match_result.num_aligned_chars)
          TreeNodePDB.objects.create(sequence_hmm = sequence_hmm,
                                      tree_node = tree_node,
                                      pdb_chain = pdb_chain)

def update_sequence_info(record,
                         seed_id,
                         seqid_of_description,
                         aligned_sequence_of_seqid,
                         num_aligned_columns,
                         sequence_of_seqid,
                         sequence_header_of_seqid,
                         assume_seed_first,
                         seed_sequence_header):
  """
  record=biopython Bio.SeqRecord.SeqRecord. 
  record.description is string (e.g., O28424_ARCFU/3-199 1,197)
  seqid_of_description - dictionary (reverse of seqs)

  """

  # Find the SEQ number associated with the record
  description = record.description
  seqid = seqid_of_description[description]

  # Find or create the aligned sequence record
  aligned_seq = record.seq.tostring()
  aligned_seguid = CheckSum.seguid(record.seq)

  aligned_sequence_of_seqid[seqid] = _sequence(AlignedSequence,
      aligned_seq, aligned_seguid)

  # Compute the number of aligned columns, if not already known
  if num_aligned_columns == 0:
    num_aligned_columns = len(aligned_seq.translate(trivial_translation,
                                                    dotlowercase))
  
  # Find or create the unaligned sequence record
  unaligned_seq = aligned_seq.translate(uppercase_translation, dotdash)
  unaligned_seguid = CheckSum.seguid(unaligned_seq)
  sequence_of_seqid[seqid] = _sequence(Sequence,
                                       unaligned_seq, unaligned_seguid)

  # Find or create the sequence_header record
  sequence_header_objects = \
      SequenceHeader.objects.filter(header__exact = description,
                                    sequence__exact = sequence_of_seqid[seqid])
  if sequence_header_objects:
    # Since the combination of header and sequence_id is constrained to be
    # unique, there can only be one
    sequence_header_of_seqid[seqid] = sequence_header_objects[0]
  else:
    # Create a new sequence_header record
    sequence_header_of_seqid[seqid] =\
                             create_sequence_header(record.id,
                                                    description,
                                                    sequence_of_seqid[seqid])
  if seed_sequence_header is None:
    if assume_seed_first:
      seed_sequence_header = sequence_header_of_seqid[seqid]
    elif seed_id is not None and seed_id == record.id:
      seed_sequence_header = sequence_header_of_seqid[seqid]
  return seed_sequence_header

def get_sequence_header_of_seqid(basename, family_accession):
  tree_with_seqids_filename = basename + ".nj.rooted.tre"
  tree_with_seqhdrs_filename = family_accession + ".nj"
  f = open(tree_with_seqids_filename)
  seqid_tree_str = f.read()
  f.close()
  seqid_tree_tokens = seqid_re.split(seqid_tree_str)
  f = open(tree_with_seqhdrs_filename)
  seqhdr_tree_str = f.read()
  f.close()
  seqhdr_tree_tokens = seqhdr_re.split(seqhdr_tree_str)
  sequence_header_of_seqid = {}
  for i in range(len(seqid_tree_tokens)):
    if i % 2 == 1:
      seqid = seqid_tree_tokens[i]
      seqhdr = seqhdr_tree_tokens[i]
      m = capturing_seqhdr_re.match(seqhdr)
      sequence_header = SequenceHeader.objects.get(id = m.group(2))
      sequence_header_of_seqid[seqid] = sequence_header
  return sequence_header_of_seqid
  

@transaction.commit_on_success
def insertMLTree(alignment_path, family_accession):
  family_id = int(family_accession[3:])
  family = Family.objects.get(id = family_id)
  basename = os.path.splitext(alignment_path)[0]
  ml_tree_filename = basename + ".fasttree.ml.rooted.tre"
  nj_tree_filename = basename + ".nj.rooted.tre"
  if os.path.exists(nj_tree_filename) and os.path.exists(ml_tree_filename):
    root = node()
    f = open(ml_tree_filename)
    treeString = f.read()
    f.close()
    treeString = treeString.translate(trivial_translation, string.whitespace)
    root.readFromTreeString(treeString, 0)
    tree = Tree.objects.create(family=family, method='ml', is_rsd_rooted=False)
    root.updateLeftId(1, 0)
    sequence_header_of_seqid = get_sequence_header_of_seqid(basename,
                                                            family_accession)
    root.createTreeNodeObjects(tree, sequence_header_of_seqid)
    other_root = family.canonical_root_node()
    other_tree_node_alignments = TreeNodeAlignment.objects.filter(
                                                      tree_node = other_root)
    for tree_node_alignment in other_tree_node_alignments:
      TreeNodeAlignment.objects.create(tree_node=root.tree_node,
                      aligned_sequence = tree_node_alignment.aligned_sequence,
                      sequence_header = tree_node_alignment.sequence_header)
    replace_seqids_by_seq_header_ids(ml_tree_filename, 
                                      family_accession + '.ml',
                                      sequence_header_of_seqid)
    return True
  else:
    return False

def insertAlignmentConservation(basename, canonical_root_node):
  aligned_column_indices = set()
  column_conserved_residue = {}
  column_score = {}
  infname = basename + '.alignmentconservation.csv'
  inf = open(infname)
  # skip header line
  for line in inf.readlines()[1:]:
    fields = line.strip().split(',')
    j = int(fields[0])
    score = float(fields[2])
    aligned_column_indices.add(j)
    column_conserved_residue[j] = fields[1]
    column_score[j] = score
  inf.close()
  update_tree_node_alignment_conservation(canonical_root_node, 
                                          aligned_column_indices,
                                          column_conserved_residue,
                                          column_score)

def insertFamilyIntoDB(alignment_path, assume_seed_first=False, seed_id=None,
    build_database_source="UniProt",
    gathering_method="FlowerPower", 
    private=False,
    family_specific_evalue_criterion=None,
    famiy_specific_sw_method=None,
    notes=None,
    build_alignment_notes_id=0,
    family_type_id='C'):

  workdir, alignment_filename = os.path.split(alignment_path)
  os.chdir(workdir)

  # Read in the alignment
  f = open(alignment_filename)
  alignments = AlignIO.parse(f, 'fasta')

  # AlignIO.parse returns a list of alignments. We only want one
  # alignment, so we take the first one (there shouldn't be any more
  # for our inputs anyhow)
  for alignment in alignments:
    break

  basename = os.path.splitext(alignment_filename)[0]
  msg = 'This file should have been created by buildFamily.py.'

  # Look for the ID mapping from SEQ num identifiers to fasta headers
  idmap_filename = basename + '.idmap'
  if not os.path.exists(idmap_filename):
    print "File %s not found. %s" % (idmap_filename, msg)
    return 1

  f = open(idmap_filename)
  idmap = cPickle.load(f)
  f.close()

  # Reverse the ID mapping
  seqid_of_description = {}
  for seqid in idmap:
    seqid_of_description[idmap[seqid]] = seqid

  # Assert assumptions from buildFamily. The idmap is a mapping between
  # both unique SEQs and unique headers
  assert(len(seqid_of_description.keys()) == len(idmap.keys()))

  sequence_of_seqid = {}
  aligned_sequence_of_seqid = {}
  sequence_header_of_seqid = {}
  seed_sequence_header = None
  num_aligned_columns = 0
  # For each sequence in the alignment, make sure records exist in the
  # sequence, sequence_header, and aligned_sequence tables
  for record in alignment:
    seed_sequence_header = update_sequence_info(record,
                         seed_id,
                         seqid_of_description,
                         aligned_sequence_of_seqid,
                         num_aligned_columns,
                         sequence_of_seqid,
                         sequence_header_of_seqid,
                         assume_seed_first,
                         seed_sequence_header)

  canonical_tree_method = ""
  root_of_method = {}
  if len(sequence_header_of_seqid) < 4:
    # There won't be an actual tree, so make a fake one
    root_of_method['trivial'] = node()
    for seqid in sequence_header_of_seqid.keys():
      child = node(seqid=seqid)
      child.branch_length = 1.0
      root_of_method['trivial'].addChild(child)
    canonical_tree_method = "trivial"
    ml_tree_filename = None
    nj_tree_filename = None
  else:
    nj_tree_filename = basename + ".nj.rooted.tre"
    ml_tree_filename = basename + ".fasttree.ml.rooted.tre"
    # TODO: We should also check that the tree files are nonempty.  If both are
    # empty we should create a trivial tree.
    if os.path.exists(ml_tree_filename):
      root_of_method['ml'] = node()
      f = open(ml_tree_filename)
      treeString = f.read()
      f.close()
      treeString = treeString.translate(trivial_translation, string.whitespace)
      root_of_method['ml'].readFromTreeString(treeString, 0)
      canonical_tree_method = "ml"
    if os.path.exists(nj_tree_filename):
      root_of_method['nj'] = node()
      f = open(nj_tree_filename)
      treeString = f.read()
      f.close()
      treeString = treeString.translate(trivial_translation, string.whitespace)
      root_of_method['nj'].readFromTreeString(treeString, 0)
      if canonical_tree_method == "":
        canonical_tree_method = "nj"

  if canonical_tree_method == "":
    print "No tree file found. %s" % msg
    return 1

  # Try to read the build date from a file, otherwise assume it is today
  build_date_filename = basename + ".build_date"
  build_date = datetime.date.today()
  if os.path.exists(build_date_filename):
    try:
      f = open(build_date_filename)
      year, month, day = [int(field) for field in f.read().strip().split('-')]
      f.close()
      build_date = datetime.date(year, month, day)
    except ValueError:
      pass

  # Try to get the build_alignment_notes, if applicable
  if build_alignment_notes_id > 0:
    build_alignment_notes = BuildAlignmentNotes.objects.get(
                                          id__exact=build_alignment_notes_id)

  # At this point we have the minimum information necessary, namely an
  # alignment and a tree, so we can go ahead and create a family

  # TODO: We should set the status to "bad" here, and then update the status to
  # "draft" at the very end.
  family = Family.objects.create(build_database_source=build_database_source,
                            build_date=build_date,
                            status="draft", # all families start out as draft
                            private=private,
                            gathering_method=gathering_method,
                            family_type_id=family_type_id,
                            partition="B", # all families start in B partition
                            )

  # Now we have created a new family accession
  family_accession = 'bpg%07d' % family.id
  print family_accession

  # Create the appropriate directories and symbolic links 
  pfacts_base_dir = '/clusterfs/ohana/bpg/pfacts'
  dir1 = os.path.join(pfacts_base_dir, family_accession[0:4])
  dir2 = os.path.join(dir1, family_accession[0:7])
  dir3 = os.path.join(dir2, family_accession)
  if not os.path.exists(dir1):
    os.mkdir(dir1)
  if not os.path.exists(dir2):
    os.mkdir(dir2)
  if not os.path.exists(dir3):
    os.chdir(dir2)
    os.symlink(workdir, family_accession)
    os.chdir(workdir)

  # Create symbolic links to files
  os.symlink(alignment_filename, family_accession + '.a2m')
  os.symlink(idmap_filename, family_accession + '.idmap')
  if ml_tree_filename is not None and os.path.exists(ml_tree_filename):
    replace_seqids_by_seq_header_ids(ml_tree_filename, 
                                      family_accession + '.ml',
                                      sequence_header_of_seqid)
  if nj_tree_filename is not None and os.path.exists(nj_tree_filename):
    replace_seqids_by_seq_header_ids(nj_tree_filename, 
                                      family_accession + '.nj',
                                      sequence_header_of_seqid)

  # Link in the build_alignment_notes
  if build_alignment_notes_id > 0:
    family.build_alignment_notes = build_alignment_notes
    family.save()

  # Create the tree objects
  tree_of_method = {}
  for method in ['trivial', 'nj', 'ml']:
    if method in root_of_method:
      tree_of_method[method] = Tree.objects.create(family=family,
                                                  method=method,
                                                  is_rsd_rooted=False
                                                  )

  canonical_tree = tree_of_method[canonical_tree_method]

  # Link the canonical tree to the family
  family.canonical_tree = canonical_tree
  family.save()

  # Link the seed sequence header to the family
  if seed_sequence_header:
    family.seed_sequence_header = seed_sequence_header
    family.save()

  # Do the modified pre-order tree traversal to find the leftIds and rightIds
  # of each of the nodes
  for method in ['trivial', 'nj', 'ml']:
    if method in root_of_method:
      root_of_method[method].updateLeftId(1, 0)

  # Create tree_node objects for each tree
  # The leaf nodes will be linked to sequence_header objects
  for method in ['trivial', 'nj', 'ml']:
    if method in root_of_method:
      root_of_method[method].createTreeNodeObjects(tree_of_method[method],
                                                  sequence_header_of_seqid)

  # Link the family alignment to the root of each tree
  # It appears redundant to link the family alignment multiple times, but this
  # is not the case.  If we subsequently run SATCHMO-JS, we will create a new
  # tree for this family with method 'satchmo-js', and the alignment linked to
  # the root of that tree will be a different alignment, namely the one output
  # by SATCHMO-JS.  We will also link the SATCHMO-JS subalignments to the
  # internal nodes of that tree, and we may make that tree the canonical tree.
  for method in ['trivial', 'nj', 'ml']:
    if method in root_of_method:
      root_node = root_of_method[method].tree_node
      for seqid in aligned_sequence_of_seqid:
        TreeNodeAlignment.objects.create(tree_node=root_node,
                            aligned_sequence=aligned_sequence_of_seqid[seqid],
                            sequence_header=sequence_header_of_seqid[seqid])

  # Now link some family data to the root of the canonical tree
  # TODO: Write a procedure for changing the canonical tree.  This procedure
  # must link all this family data to the root of the new canonical tree.
  canonical_root_node = root_of_method[canonical_tree_method].tree_node

  # Link the family HMMs
  sam_hmm_filename = basename + '.mod'
  if os.path.exists(sam_hmm_filename):
    os.symlink(sam_hmm_filename, family_accession + '.mod')
    sam_hmm = HMM.objects.create(length=num_aligned_columns,
                                hmm_type = 'SAM',
                                method = 'w0.5',
                                tree_node = canonical_root_node)


  hmmer_hmm_filename = basename + '.hmm'
  if os.path.exists(hmmer_hmm_filename):
    # Rewrite the HMMER hmm so its name is the family id
    inf = open(hmmer_hmm_filename)
    hmm_lines = inf.readlines()
    inf.close()
    outf = open((family_accession + '.hmm'), "w")
    for line in hmm_lines:
      if line[0:4] == 'NAME':
        outf.write("NAME  %s\n" % family_accession)
      else:
        outf.write(line)
    outf.close()
    hmmer_hmm = HMM.objects.create(length=num_aligned_columns,
                                hmm_type = 'HMMER3',
                                method = 'hmmbuild',
                                tree_node = canonical_root_node)

  # Link the family consensus sequence
  consensus_sequence_filename = basename + '.con.fa'
  if os.path.exists(consensus_sequence_filename):
    os.symlink(consensus_sequence_filename, family_accession + '.con.fa')
    f = open(consensus_sequence_filename)
    record = list(SeqIO.parse(f, "fasta"))[0]
    f.close()
    consensus_seguid = CheckSum.seguid(record.seq)
    # It's not inconceivable that the consensus sequence is already in the
    # sequence table.  Find the sequence record for this consensus sequence, or
    # create it if it isn't there

    consensus_sequence = _sequence(Sequence, 
                                   record.seq.tostring(), consensus_seguid)

    # Link the consensus sequence to the family hmm
    hmmer_hmm_consensus = HMM_Consensus.objects.create(hmm = hmmer_hmm,
                                                  sequence = consensus_sequence)
    # Link the consensus sequence to the canonical tree root.
    # It appears redundant that the sequence record is linked here, since the
    # sequence record is already linked to the hmm_consensus record.  The
    # reason for linking the sequence record directly to the
    # tree_node_consensus record is that the consensus sequence might not come
    # from an HMM; it might be derived directly from the alignment instead.  In
    # that case there would have been no hmm_consensus record linked to the
    # tree_node_consensus record.  But there must always be a sequence record
    # linked to every tree_node_consensus record.
    canonical_root_consensus = TreeNodeConsensus.objects.create(
                                tree_node = canonical_root_node,
                                sequence = consensus_sequence,
                                method = 'hmm',
                                hmm_consensus = hmmer_hmm_consensus)

  # Insert the alignment conservation
  os.symlink(basename + '.alignmentconservation.csv', 
              family_accession + '.alignmentconservation.csv')
  insertAlignmentConservation(family_accession, canonical_root_node)

  # Link the PFAM domains
  insertPFAMPredictionsIntoDB(consensus_sequence, basename)

  # Link the signal peptide and transmembrane prediction
  phobius_filename = basename + '.phobius'
  if os.path.exists(phobius_filename):
    insertPhobiusPredictionsIntoDB(canonical_root_node, phobius_filename)

  # Link the homologous PDB structures
  insertPDBPredictionsIntoDB(hmmer_hmm, canonical_root_node, basename)

  if os.path.exists(build_date_filename):
    os.symlink(build_date_filename, family_accession + '.build_date')

def main():
  parser = OptionParser(usage='%prog [options] alignment_path seed_id')
  parser.add_option('--assume_seed_first', action='store_true',
                dest='assume_seed_first', default=False,
                help="Assume the first sequence in the alignment is the seed.")
  parser.add_option('--no_assume_seed_first', action='store_false',
          dest='assume_seed_first', default=False,
          help="Don't assume the first sequence in the alignment is the seed.")
  parser.add_option('--force_no_seed_id', action="store_true",
                    dest='force_no_seed_id', default=False,
                    help="Insert family without seed sequence identifier")
  parser.add_option('--build_database_source', dest='build_database_source',
                    default="UniProt",
                    help="Database from which homologs were gathered")
  parser.add_option('--private', action='store_true', dest='private',
                    default=False,
                    help="Make the family private")
  parser.add_option('--public', action='store_false', dest='private',
                    default=False,
                    help="Make the family public")
  parser.add_option('--gathering_method', dest='gathering_method',
                    default='FlowerPower',
                    help="Method used to gather homologs")
  parser.add_option('--notes', dest='notes', default='', 
                    help="Family-specific notes")
  parser.add_option('--build_alignment_notes_id',
                    dest='build_alignment_notes_id',
                    default=0,
                    help="Id of the build aligment notes entry in the database")
  parser.add_option('--family_type',
                    dest='family_type',
                    default='C',
                    help="Family Type (usually either 'C' or 'G')")
  (options, args) = parser.parse_args()

  if len(args) < 1:
    parser.error('Must supply the path to an alignment file.')

  if not options.force_no_seed_id and len(args) < 2:
    parser.error('Must supply the identifier of the seed sequence '
                + '(use --force_no_seed_id to override)')

  alignment_path = args[0]
  if options.force_no_seed_id:
    seed_id = None
  else:
    seed_id = args[1]
  status_code = insertFamilyIntoDB(alignment_path,
            assume_seed_first=options.assume_seed_first,
            seed_id = seed_id,
            build_database_source=options.build_database_source,
            gathering_method=options.gathering_method,
            private=options.private,
            notes=options.notes,
            build_alignment_notes_id=options.build_alignment_notes_id,
            family_type_id=options.family_type,
            )
  sys.exit(status_code)

if __name__ == '__main__':
  main()
