#!/usr/bin/env python

import os, sys

from pfacts003.phylofacts.models import Tree
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

  trees = Tree.objects.filter(family__id = family_id)
  had_error = False
  for tree in trees:
    filename = '%s_subst.%s' % (family_accession, tree.method)
    try:
      outf = open(filename, "w")
      tree.write_newick_to_handle(outf)
      outf.close()
      sys.stdout.write("Wrote %s\n" % filename)
    except IOError:
      sys.stderr.write("IOError while attempting to write %s\n" % filename)
      had_error = True

  if had_error:
    sys.exit(1)
  sys.exit(0)

if __name__ == '__main__':
  main()
