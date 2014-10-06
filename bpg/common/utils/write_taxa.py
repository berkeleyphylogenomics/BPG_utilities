#!/usr/bin/env python

import os, sys

from pfacts003.phylofacts.models import Family, TreeNodeAlignment, SequenceHeader
from bpg.common.utils.dir_of_family import get_dir_of_family_accession

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <family_accession>" % sys.argv[0]
    sys.exit(0)

  family_accession = sys.argv[1]
  try:
    family_id = int(family_accession[4:])
  except ValueError:
    print "Usage: %s <family_accession>" % sys.argv[0]
    sys.exit(1)

  family_dir = get_dir_of_family_accession(family_accession)
  if not os.path.exists(family_dir):
    print "Family %s directory not found" % family_accession
    sys.exit(1)
  os.chdir(family_dir)

  family = Family.objects.get(id = family_id)
  tree_node_alignment_objects = TreeNodeAlignment.objects.filter(
                                  tree_node = family.canonical_root_node())
  for obj in tree_node_alignment_objects:
    if obj.sequence_header.uniprot is not None:
      print "%s\t%d" % (obj.sequence_header.get_escaped_header_for_tree(),
                        obj.sequence_header.uniprot.taxon.id)

if __name__ == '__main__':
  main()
