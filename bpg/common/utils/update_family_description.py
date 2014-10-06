#!/usr/bin/env python

from pfacts003.phylofacts.models import Family, TreeNodeName, OrthologTypes
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if len(args[0]) != 10 or args[0][0:3] != 'bpg':
    opt_parser.error('Argument must be a bpg accession like bpg0164198')
  bpg_accession = args[0]
  try:
    family_id = int(bpg_accession[3:])
  except ValueError:
    opt_parser.error('Argument must be a bpg accession like bpg0164198')
  family = Family.objects.get(id = family_id)
  root = family.canonical_root_node()
  if not root.canonical_tree_node_name:
    name = root.get_description(ortholog_type = OrthologTypes.PHOG_T_Custom,
                                threshold = 10000.0)
    tree_node_name = TreeNodeName.objects.create(tree_node = root,
                                          name = name,
                                          source = 'uniprot',
                                          method = 'plurality SwissProt 10')
    root.canonical_tree_node_name = tree_node_name
    root.save()

if __name__ == '__main__':
  main()

