from __future__ import division
from files import *
import cPickle
import logging
import os
import re
import shlex
import StringIO
import subprocess
import time
import sys
from utils import remove_duplicates

from Bio import SeqIO, Blast
from Bio.SeqUtils.CheckSum import seguid
from Bio.Blast import NCBIXML
from Bio.SubsMat import MatrixInfo
from django.utils.safestring import mark_safe
from django.db.models import Q

from pfacts003.utils.id_patterns import fasta_re
from bpg.common.parsers.hmm_parsers import \
    parse_results_of_hmmsearch_or_hmmscan
from pfacts003.phylofacts.models import TreeNode, SequenceHeader, UniProt, \
    Family, HMM, getTreeNodeQuerySetForTreeNodes
from pfacts003.utils.id_patterns import *
from pfacts003.utils.messages import unrecognized_identifier, not_in_phog
from pfacts003.utils.extract_sequence_from_fasta \
  import extract_sequence_from_fasta
from pfacts003.phylofacts.models import OrthologTypes
from pfacts003.utils.links import make_pfam_links, alignment_hmm_pairwise_link

logging.basicConfig(level = logging.INFO)

def getTimeStr(starttime, endtime):

    """Returns time interval formatted in minutes and seconds

    This function formats a time interval into minutes and seconds, so we can
display the elapsed time for users
"""
    elapsed_secs = endtime - starttime

    if elapsed_secs > 60:
        mins = int(elapsed_secs / 60)
        elapsed_secs = elapsed_secs - (mins * 60)

        if mins == 1:
            return '[1 min, %.2f secs]' % elapsed_secs

        return '[%d mins, %.2f secs]' % (mins, elapsed_secs)

    return '[%.2f secs]' % elapsed_secs

def fasta_file_to_record(fasta_file_name):
  f = open(fasta_file_name, "rU")
  records = list(SeqIO.parse(f, "fasta"))
  return(records[0])

ghg_family_hmm_file = '/clusterfs/ohana/bpg/pfacts/phylofacts3_GHG_090811.hmms'
pfam_family_hmm_file = '/clusterfs/ohana/bpg/pfacts/phylofacts3_domain_090811.hmms'
def approximately_matching_families_from_fasta(record, context=None,
                                                extra_context=None,
                                                e_value_cutoff = 1e-10,
                                                pc_identity_cutoff = 70.0):
  querylen = len(record.seq.tostring())
  (hmmscan_input, base_name) = create_file(PHOG_BLAST_DIR, "_input.fa")
  f = open(hmmscan_input, "w")
  SeqIO.write([record], f, "fasta")
  f.close()
  
  def run_and_parse_hmmscan(family_type_str, family_hmm_file):
    hmmscan_results_filename = '_'.join([base_name, family_type_str, "hmmscan.out"])
    hmmscan_results_table_filename = '_'.join([base_name, family_type_str, "hmmscan.domtblout"])
    hmmscan_output = PHOG_BLAST_DIR + hmmscan_results_filename
    hmmscan_tbloutput = PHOG_BLAST_DIR + hmmscan_results_table_filename    
    starttime = time.time()
    logging.info('Running %s hmmscan' % family_type_str)
    sys.stdout.flush()
    p = subprocess.Popen(['/clusterfs/ohana/software/bin/hmmscan',
                        '-o', hmmscan_output,
                        '-E', str(e_value_cutoff),
                        '--incE', str(e_value_cutoff),
                        '--domtblout', hmmscan_tbloutput,
                        family_hmm_file, 
                        hmmscan_input])
    output, error = p.communicate(input=None)
    endtime = time.time()
    logging.info('Done. %s' % getTimeStr(starttime, endtime))
    starttime = time.time()
    logging.info('Parsing results')
    sys.stdout.flush()
    # For GHG families, require that the hit align globally to the query
    if family_type_str == "ghg":
        min_sequence_length = int(.7*len(record.seq.tostring()))
    else:
        min_sequence_length = 10
    hmmscan_results_by_query = parse_results_of_hmmsearch_or_hmmscan.parse(
                        PHOG_BLAST_DIR + hmmscan_results_filename,
                        e_value_cutoff, min_sequence_length, 10)
    endtime = time.time()
    logging.info('Done. %s' % getTimeStr(starttime, endtime))
    starttime = time.time()
    logging.info('Parsing domtbl')
    sys.stdout.flush()
    ordered_families = []
    included_families = set()
    f = open(PHOG_BLAST_DIR + hmmscan_results_table_filename)
    for line in f.readlines():
      if len(line) > 0 and line[0] != '#':
        family = line.split()[0]
        if family not in included_families:
          ordered_families = ordered_families + [family]
          included_families.add(family)
    f.close()
    endtime = time.time()
    logging.info('Done. %s' % getTimeStr(starttime, endtime))
    for query in hmmscan_results_by_query.hit_result_of_name_of_query:
      hmmscan_results = hmmscan_results_by_query.hit_result_of_name_of_query[
                                                                        query]
      break
    highest_scoring_match_number = {}
    best_bit_score = {}
    best_match = {}
    starttime = time.time()
    logging.info('Finding best matches')
    sys.stdout.flush()
    for family_name in ordered_families:
      if family_name not in hmmscan_results:
        continue
      family = Family.objects.get(id__exact = int(family_name[3:]))
      if family.is_global_homology():
        family_type = 1
      else:
        family_type = 0
      highest_scoring_match_number[family.id] = -1
      best_bit_score[family.id] = -10000.0
      for match_number in hmmscan_results[family_name].matches:
        match_result = hmmscan_results[family_name].matches[match_number]
        if match_result.bit_score > best_bit_score[family.id]:
          best_bit_score[family.id] = match_result.bit_score
          highest_scoring_match_number[family.id] = match_number
      if highest_scoring_match_number[family.id] >= 0:
        best_match[family.id] \
          = hmmscan_results[family_name].matches[
                                        highest_scoring_match_number[family.id]]
      root = family.canonical_root_node()
      hmm = HMM.objects.get(tree_node = root, hmm_type = 'HMMER3')
    endtime = time.time()
    logging.info('Done. %s' % getTimeStr(starttime, endtime))
    output_file = os.path.join(PHOG_BLAST_DIR, base_name + '_%s_output.pkl' % family_type_str)
    f = open(output_file, 'w') 
    cPickle.dump(best_match, f)
    f.close()
    return (ordered_families, best_match)
  (ordered_ghg_families, best_ghg_match) = run_and_parse_hmmscan('ghg', ghg_family_hmm_file)
  (ordered_pfam_families, best_pfam_match) = run_and_parse_hmmscan('pfam', pfam_family_hmm_file)
  if not extra_context is None:
    extra_context['base_name'] = base_name
  return( ( (ordered_ghg_families, best_ghg_match), (ordered_pfam_families, best_pfam_match) ) )


def get_family_data_for_sequence(defline, sequence):
    """
    Calculates and returns the data needed for sequence search. Also see 
    def sequence_search in phylofacts/views.py.
    """
    context = {}
    extra_context = {}
    record = get_seq_record_from_defline_and_sequence(defline, sequence)
    if record is None:
       return {}
    ((ghg_family_ids, best_ghg_match), (pfam_family_ids, best_pfam_match), error) = get_approximately_matching_family_query_set(record, context, extra_context)
    context.update(extra_context)
    # Removed 'and (pfam_family_ids == [] or best_pfam_match == {}) from the code below, since we are only
    # displaying GHG books in this iteration. Keeping it would mean that if there are no GHG matches in
    # the family, but PFAM matches are available, then the 'No families' message will not be printed.
    if (ghg_family_ids == [] or best_ghg_match == {}) or (error != None):
        return {}
    def fetch_ordered_families_from_ids(family_ids):

      family_ids = family_ids[:20]

      # Make a single round trip to the database and ...
      families = Family.objects.filter(id__in=family_ids)
      family_ids_dict = {}

      # ... re-order to revert to the original order of family_ids.
      for family in families:
          family_ids_dict[family.id] = family
      families = [family_ids_dict[family_id] for family_id in family_ids] 
      return families

    all_ghg_families = fetch_ordered_families_from_ids(ghg_family_ids)
    ghg_families = [family for family in all_ghg_families 
    if family.id in best_ghg_match and 
    (best_ghg_match[family.id].seq_to - 
      best_ghg_match[family.id].seq_from) 
      >= 0.7 * family.get_alignment_length()]
    # If the above alignment size filter does not produce any results, then return 'Not found' message
    if not ghg_families: return {}
    pfam_families = fetch_ordered_families_from_ids(pfam_family_ids)
    # Populate the context dict...
    context['sequence'] = sequence
    context['query_title'] = record.description
    def best_match_context(families, best_match):
      return [{'description':family.canonical_root_node().get_description(ortholog_type=OrthologTypes.PHOG_T_Custom, threshold=10000.0),
        'family_link':mark_safe(family.get_link()),
        'query_title':record.description,
        'alignment_link':alignment_hmm_pairwise_link('Alignment',
             family.id, record.description,
             len(record.seq.tostring()),
             extra_context['base_name']),
        'best_match':best_match[family.id]} for family in families if family.id in best_match]
    context['best_ghg_matches'] = best_match_context(ghg_families, best_ghg_match)
    context['best_pfam_matches'] = best_match_context(pfam_families, best_pfam_match)
    return context

def approximately_matching_seq_headers_from_fasta(record,
                                                  context = None,
                                                  extra_context = None,
                                                  e_value_cutoff = 1,
                                                  pc_identity_cutoff = 70.0):

  def parse_blast_results(_filter_by_identity=False):
    seqhdr_alignment_dict = {}
    f = open(blast_output, "r")
    blast_records = NCBIXML.parse(f)
    
    for blast_record in blast_records:
      for alignment in blast_record.alignments:
        m = re.search("bpgprotseq\d+", alignment.title)
        if m:
          seqhdr_id = str(int(m.group(0)[10:]))
          best_hsp = None
          best_e_value = 1000.0
  
          for hsp in alignment.hsps:
            if hsp.expect <= e_value_cutoff and hsp.expect < best_e_value:
              best_e_value = hsp.expect
              best_hsp = hsp
              
          if best_hsp:
            percent_id = (best_hsp.identities/len(best_hsp.query))*100.0
            if percent_id >= pc_identity_cutoff or not _filter_by_identity:
              seqhdr_alignment_dict[seqhdr_id] \
                  = {'percent_identity': percent_id,
                     'score': best_hsp.score,
                     'e_value': best_hsp.expect}
    f.close()
    return seqhdr_alignment_dict                
  
  skip_blast = record is None
  if skip_blast: 
    record = fasta_file_to_record("/home/meacham/fa/test.fa")
  qlen = len(str(record.seq))
  (blast_input, base_name) = create_file(PHOG_BLAST_DIR, "_input.fa")
  blast_results_file_name = base_name + "_results.xml"
  blast_output = PHOG_BLAST_DIR + blast_results_file_name
  f = open(blast_input, "w")
  SeqIO.write([record], f, "fasta")
  f.close()
  
  blast_db  = BLAST_DB_DIR + 'bpgseqs_in_phogs/bpgseqs_in_phogs'
  if not skip_blast:
    os.system("/clusterfs/ohana/software/bin/blastall -p blastp -m 7 -F F -d %s -i %s -o %s" 
            % (blast_db, blast_input, blast_output))

  seqhdr_alignment_dict = parse_blast_results(True)
  
  if len(seqhdr_alignment_dict) != 0:
    if not context is None:
      context['blast_results_message'] = mark_safe("BLAST matches to sequences in "
        + "PhyloFacts Orthology Groups")
      context['blast_results_explanation'] = mark_safe("Please note:&nbsp; "
        + "Sequences shown below are homologous, but may not be orthologous. "
        + "If an exact match exists (100% ID), <span style='color: red'>click "
        + "on the orthologs link</span> to find orthologs to that "
        + "sequence based on phylogenetic analysis. Caution should be exercised "
        + "for matches with less than 80% ID as these may be paralogous.")
  else:

    seqhdr_alignment_dict = parse_blast_results(False)
  
    if len(seqhdr_alignment_dict) != 0:
      if not context is None:
        context['blast_results_message'] = \
          ('No sequences with better than %.f%% identity were found in any PhyloFacts orthology group. ' \
          + 'These sequences were found with a BLAST E-value of better than %.1e. However, we warn ' \
          + 'that orthology inferences based on sequences with only moderate ' \
          + 'sequence identity and/or fractional overlap can be misleading.') \
          % (pc_identity_cutoff, e_value_cutoff)
    elif not context is None:
      context['blast_results_message'] = \
        ('No sequences with an E-value of better than %.1e were found in ' 
         + 'any PhyloFacts orthology group.') % e_value_cutoff

  if not extra_context is None:
    extra_context['blast_results_file_name'] = blast_results_file_name
  return(seqhdr_alignment_dict)
       
def get_approximately_matching_family_query_set(record, context = None, extra_context = None):
  ((ordered_ghg_families, best_ghg_match), (ordered_pfam_families, best_pfam_match)) \
    = approximately_matching_families_from_fasta(record, context, extra_context)
  if len(ordered_ghg_families) == 0 and len(ordered_pfam_families) == 0:
    return (([], {}), ([], {}), None)
  ghg_family_ids = [int(family_accession[4:]) for family_accession in
                ordered_ghg_families]
  pfam_family_ids = [int(family_accession[4:]) for family_accession in
                ordered_pfam_families]
  return ((ghg_family_ids, best_ghg_match), (pfam_family_ids, best_pfam_match), None)

def get_approximately_matching_sequence_query_set(record, context = None, extra_context = None,
                                    orderByTaxonomicDistribution = False, number_to_retain=20):
  (seqhdr_alignment_dict) \
      = approximately_matching_seq_headers_from_fasta(record, context, extra_context)   
      
  seqhdr_ids = []
  for seqhdr_id in seqhdr_alignment_dict:
    seqhdr_ids.append(seqhdr_id) # May contain duplicate seqhdr_ids

  approximate_matches = TreeNode.objects.filter(
    sequence_header__id__in=seqhdr_ids)\
      .exclude(sequence_header__uniprot__in=[0])\
      .exclude(sequence_header__uniprot__accession__exact=None)

  if not bool(approximate_matches):
    return (None, None, not_in_phog('Sequence'))

# Sort results by BLAST score (descending, hence negatives) 
  matches = [];
  for tree_node in approximate_matches:
    seqhdr_id = str(tree_node.sequence_header.id);
    matches.append((-(seqhdr_alignment_dict[seqhdr_id]['score']), str(tree_node.sequence_header), tree_node))
  matches.sort()

# Keep only one tree_node (corresponding to best score) for each sequence_header
  tree_node_list = []
  included = []
  for item in matches:
    if not item[1] in included:
      tree_node_list.append(item[2])
      included.append(item[1])
    
  tree_node_list = tree_node_list[0:number_to_retain] # Retain only the best scores

# Obtain new query set
  approximate_matches = getTreeNodeQuerySetForTreeNodes(tree_node_list)

  approximate_match_dict = {}
  for approximate_match in approximate_matches:
    id = approximate_match.sequence_header.id
    if id not in approximate_match_dict:
      approximate_match_dict[id] = set()
    if approximate_match.tree.family.is_global_homology():
      family_type = 1
    else:
      family_type = 0
    num_orthologs \
      = approximate_match.tree.family.canonical_root_node().right_id
    approximate_match_dict[id].add((family_type, num_orthologs, approximate_match))
  best_approximate_matches = set()
  for id in approximate_match_dict.keys():
    candidate_ids = list(approximate_match_dict[id])
    candidate_ids.sort(reverse=True)
    family_type, num_orthologs, approximate_match = candidate_ids[0]
    best_approximate_matches.add(approximate_match)
  if orderByTaxonomicDistribution: # ordering may be superceded by reordering in the datagrid
    approximate_matches = getTreeNodeQuerySetForTreeNodes(
                          best_approximate_matches).distinct().order_by(
                          'sequence_header__uniprot__taxon__left_id')
  else:
    approximate_matches = getTreeNodeQuerySetForTreeNodes(
                        best_approximate_matches).distinct().order_by(
                      'sequence_header__uniprot__taxon__scientific_name')
  return (approximate_matches, seqhdr_alignment_dict, None)

def get_all_uniprot_objects_for_seq_record(record):
    aa_seguid = seguid(record.seq)
    return UniProt.objects.filter(seguid=aa_seguid).all()

def get_single_uniprot_object_from_seq_record(record):
    objects = get_all_uniprot_objects_for_seq_record(record)
    if len(objects) == 0:
        return None

    # TODO: replace these 2 with function calls...
    # We really want to match by uniprot_identifier if possible.
    for object in objects:
        uniprot_identifier = object.uniprot_identifier
        if uniprot_identifier and (uniprot_identifier in record.id):
            return object.uniprot_identifier

    # We consider taxons only if there is no identifier match.
    for object in objects:
        uniprot_taxon = object.taxon
        if (uniprot_taxon and getattr(uniprot_taxon, 'scientific_name', False) and 
            uniprot_taxon.scientific_name in record.description):
            return (object.uniprot_identifier, record)

    # Return the first record if nothing else succeeds.
    return objects[0].uniprot_identifier
 
def get_seq_record_from_defline_and_sequence(defline, sequence):
    proper_seq = '\n'.join([defline, sequence.encode('ascii')])
    handle = StringIO.StringIO(proper_seq)
    record =  SeqIO.parse(handle, "fasta").next()
    return record


# TODO: deprecated; this is not used in Phylofacts; it is
# used in PHOG and when it is replaced in PHOG by some of the 
# other defs defined here, this should be removed.
# - ST, 2/14/11
def get_uniprot_id_from_fasta(fasta_seq):
  """Retrieve UniProt identifier from a FASTA sequence
  
  Returns -- tuple (uniprot identifier, record)
  
  """
  # TODO: assume properly formatted sequence; don't do validation here
  (defline, sequence, errors) = extract_sequence_from_fasta(fasta_seq)
  if errors:
      return (None, None)
  proper_seq = '\n'.join([defline, sequence.encode('ascii')])
  handle = StringIO.StringIO(proper_seq)
  record = SeqIO.parse(handle, "fasta").next()
  record_id = record.id
  header = record.description
  aa_seguid = seguid(record.seq)

  objects = UniProt.objects.filter(seguid=aa_seguid).all()

  if len(objects) == 0:
      return (None, record)

  for object in objects:
      uniprot_identifier = object.uniprot_identifier
      # Is this a holdover from the chained headers used previously?
      if uniprot_identifier and (uniprot_identifier in record_id):
          return (object.uniprot_identifier, record)

  for object in objects:
      uniprot_taxon = object.taxon
      if uniprot_taxon and getattr(uniprot_taxon, 'scientific_name', False) and uniprot_taxon.scientific_name in header:
          return (object.uniprot_identifier, record)

  return (objects[0].uniprot_identifier, record)
  
def get_alignment_blast_pairwise(context):
  f = open_location_file(context['location'], context['file'])
  if not f is None:
    blast_records = NCBIXML.parse(f)
    for blast_record in blast_records:
      query_title = blast_record.query
      query_length = blast_record.query_length
      for alignment in blast_record.alignments:
        subject_title = alignment.title
        subject_length = alignment.length
        m = re.search("bpgprotseq\d+", subject_title)
        if m:
          this_identifier = m.group(0)
          if this_identifier == context['identifier']:
            protseq_id = this_identifier[10:]
            best_score = -1000.0
            for hsp in alignment.hsps:
              if hsp.score > best_score:
                best_score = hsp.score
                best_hsp = hsp
            identities = best_hsp.identities
            positives = best_hsp.positives
            aligned_length = len(best_hsp.query)
            blosum62 = MatrixInfo.blosum62
            best_consensus_seq = ''
            for i in range(len(best_hsp.query)):
              if best_hsp.query[i] == best_hsp.sbjct[i]:
                best_consensus_seq += best_hsp.query[i]
              else:
                try:
                  positive_score = blosum62[(best_hsp.query[i], best_hsp.sbjct[i])] > 0
                except KeyError, e:
                  try:
                    positive_score = blosum62[(best_hsp.sbjct[i], best_hsp.query[i])] > 0
                  except KeyError, e:
                    positive_score = False
                if positive_score:
                  best_consensus_seq += '+'
                else:
                  best_consensus_seq += ' '
            percent_id = '%.1f' % ((identities/aligned_length) * 100.0)
            percent_pos = '%.1f' % ((positives/aligned_length) * 100.0)
            
            if not query_length is None:
              context['query_length'] = mark_safe(query_length)
            if not subject_length is None:
              context['subject_length'] = mark_safe(subject_length)
            context['query_region'] = mark_safe('%d-%d' % (best_hsp.query_start, best_hsp.query_start+aligned_length-1))
            context['subject_region'] = mark_safe('%d-%d' % (best_hsp.sbjct_start, best_hsp.sbjct_start+aligned_length-1))
            context['identities'] = mark_safe('(%d/%d) %s%%' % (identities, aligned_length, percent_id))
            context['positives'] = mark_safe('(%d/%d) %s%%' % (positives, aligned_length, percent_pos))
            context['score'] = '%d' % best_hsp.score
            context['e_value'] = '%.1e' % best_hsp.expect
            context['query_seq'] = mark_safe(best_hsp.query.replace(' ', '&nbsp;'))
            context['consensus_seq'] = mark_safe(best_consensus_seq.replace(' ', '&nbsp;'))
            context['subject_seq'] = mark_safe(best_hsp.sbjct.replace(' ', '&nbsp;'))
            context['query_title'] = query_title
            headers = SequenceHeader.objects.filter(id__exact=protseq_id)
            if headers:
              context['subject_title'] = headers[0].identifier()
              context['subject_description'] = headers[0].description()
              context['subject_species'] = str(headers[0].uniprot.taxon)
            return

    context['error_message'] = 'The sequence identifier "%s" was not found in BLAST file "%s"' \
      % (context['identifier'], context['file'])
    f.close()
  else:
    context['error_message'] = 'The BLAST file "%s" was not found' % context['file']
  return
