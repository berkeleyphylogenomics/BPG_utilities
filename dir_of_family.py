#!/usr/bin/env python

import sys

def get_dir_of_family_accession(bpg_accession):
  return '/clusterfs/ohana/bpg/pfacts/%s/%s/%s' % (bpg_accession[0:4],
                                                  bpg_accession[0:7],
                                                  bpg_accession)

def get_dir_of_family_id(family_id):
  return get_dir_of_family_accession('bpg%07d' % family_id)

def main(): 
  if len(sys.argv) < 2:
     print "Usage: %s bpg_accession" % sys.argv[0]
     sys.exit(0)

  print get_dir_of_family_accession(sys.argv[1])

if __name__ == '__main__':
  main()
