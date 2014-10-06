#!/usr/bin/env python

import os, re, sys, glob
from pfacts003.phylofacts.models import *
from bpg.common.utils.dir_of_family import get_dir_of_family_accession
from Bio import SeqIO
from Bio.SeqUtils import CheckSum
from pfacts003.utils.id_patterns import uniprot_accession_re1, \
  uniprot_accession_re2, gi_re, swissprot_desc_pat

def usage():
  print "Usage: %s <family_accession>" % sys.argv[0]

def main():
  if len(sys.argv) < 2:
    usage()
    sys.exit(0)

  family_accession = sys.argv[1]
  try:
    family_id = int(family_accession[3:])
  except ValueError:
    usage()
    sys.exit(1)

  try:
    family = Family.objects.get(id = family_id)
    if family.status == "bad":
      raise Family.DoesNotExist
  except Family.DoesNotExist:
    print "No family found with accession %s" % family_accession
    sys.exit(1)

  family_dir = get_dir_of_family_accession(family_accession)
  seed_path = os.path.join(family_dir, "seed.fa")
  if not os.path.exists(seed_path):
    if os.path.realpath(family_dir).find('TreeFam') >= 0:
      os.chdir(family_dir)
      possible_seed_files = glob.glob("*_HUMAN*.fa")
      candidates = set()
      swissprot_desc_re = re.compile('^%s$' % swissprot_desc_pat)
      for file in possible_seed_files:
        basename = os.path.splitext(file)[0]
        components = basename.split('_')
        if len(components) < 2 or components[1] != 'HUMAN':
          continue
        if swissprot_desc_re.match(components[0]) is None and \
            uniprot_accession_re1.match(components[0]) is None and \
            uniprot_accession_re2.match(components[0]) is None:
          continue
        if len(components) > 2:
          if len(components) != 4:
            continue
          try:
            start = int(components[2])
          except ValueError:
            continue
          try:
            end = int(components[3])
          except ValueError:
            continue
        candidates.add(file)
      if len(candidates) != 1:
        print "Seed file for family %s missing" % family_accession
        sys.exit(1)
      seed_path = os.path.join(family_dir, list(candidates)[0])
    else:
      print "Seed file for family %s missing" % family_accession
      sys.exit(1)

  f = open(seed_path)
  seed_record = SeqIO.parse(f, "fasta").next()
  f.close()
  seed_seguid = CheckSum.seguid(seed_record.seq)

  seed_id = seed_record.id.strip('lcl|')
  print "%s: FlowerPower seed id %s" % (family_accession, seed_id)
  seed_accession = None
  recognizing_regexp = None
  # uniprot_accession_re1 recognizes a UniProt accession only if it is the
  # whole string, not if it is a substring
  for regexp in [re.compile(uniprot_accession_pat1),
                  re.compile(uniprot_accession_pat2), gi_re]:
    m = regexp.search(seed_id)
    if m:
      seed_accession = m.group()
      recognizing_regexp = regexp
      break
    
  if seed_accession is None:
    print "Could not parse accession from seed id"
    sys.exit(1)

  sequences = Sequence.objects.filter(seguid = seed_seguid)
  sequence_headers = SequenceHeader.objects.filter(sequence__in = sequences)
  possible_sequence_headers = set()
  for sequence_header in sequence_headers:
    m = recognizing_regexp.search(sequence_header.header)
    if m:
      accession = m.group()
      if accession == seed_accession:
        if len(sequence_header.header) >= 4 and \
            sequence_header.header[0:4] == 'lcl|':
          possible_sequence_headers = set([sequence_header])
          break
        if sequence_header.header.find('|') < 0:
          possible_sequence_headers = set([sequence_header])
          break
        possible_sequence_headers.add(sequence_header)
  if len(possible_sequence_headers) > 1:
    alns = TreeNodeAlignment.objects.filter(tree_node =
            family.canonical_root_node(), 
            sequence_header__in = possible_sequence_headers)
    possible_sequence_headers = set([aln.sequence_header for aln in alns])
    
  print "%s: Found %d possible sequence headers" % (family_accession,
                                                len(possible_sequence_headers))
  for seqhdr in possible_sequence_headers:
    print "%s: possible sequence header %s" % (family_accession, seqhdr.header)

  if len(possible_sequence_headers) == 1:
    seed_sequence_header = list(possible_sequence_headers)[0]
    print "Assigning seed sequence header id %d to family %s" \
        % (seed_sequence_header.id, family_accession)
    family.seed_sequence_header = seed_sequence_header
    family.save()

if __name__ == '__main__':
  main()
