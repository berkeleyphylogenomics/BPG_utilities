#!/usr/bin/env python

import sys
from pfacts003.utils.id_patterns import *
from pfacts003.phylofacts.models import *
from Bio import SeqIO
from Bio.SeqUtils import CheckSum
import urllib2

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <uniprot_accession>" % sys.argv[0]
    sys.exit(0)
  uniprot_accession = sys.argv[1]
  if not uniprot_accession_re1.match(uniprot_accession) and not \
      uniprot_accession_re2.match(uniprot_accession):
    print "The argument must be a valid UniProt accession"
    sys.exit(1)
  try:
    response = urllib2.urlopen('http://www.uniprot.org/uniprot/%s.fasta' 
                                % uniprot_accession)
  except urllib2.HTTPError:
    print "Unable to download sequence from UniProt"
    sys.exit(1)

  record = SeqIO.parse(response, 'fasta').next()
  seguid = CheckSum.seguid(record.seq)
  sequence_objects = Sequence.objects.filter(seguid__exact = seguid)
  if sequence_objects:
    tree_node_alignment_objects = TreeNodeAlignment.objects.filter(
                              sequence_header__sequence__in = sequence_objects)
    if tree_node_alignment_objects:
      families = set([obj.tree_node.tree.family for obj in
                    tree_node_alignment_objects])
      for family in families:
        print family.get_accession()
    else:
      print "There are no families containing this sequence."
  else:
    print "This sequence is not in the PhyloFacts 3 database"
  
if __name__ == '__main__':
  main()
