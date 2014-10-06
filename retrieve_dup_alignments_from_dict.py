#!/usr/bin/python

import os, sys, cPickle
from matchmaker.shmm_shmm_lib import *

def main():
  if len(sys.argv) < 4:
    print "Usage: %s <seedX_id> <seedY_id> <alignment_filename>" % sys.argv[0]
    print "Using the output of "
    print "  print_nr_alignments_with_dict.py <seedX_id> <seedY_id>"
    print "prints all alignments  of <seedX_id> with <seedY_id> which are " \
          + "identical to the one in <alignment_filename>."
    sys.exit(0)
  seedX_id = sys.argv[1]
  seedY_id = sys.argv[2]
  alignment_filename = sys.argv[3]
  pklfp = open("%s_%s_alignment_dict.pkl" % (seedX_id, seedY_id))
  (dups_of_alignment, nr_alignments) = cPickle.load(pklfp)
  pklfp.close()
  alignment_filenames = dups_of_alignment[alignment_filename]
  for dup_alignment_filename in alignment_filenames:
    print dup_alignment_filename
        
if __name__ == '__main__':
  main()

