#!/usr/bin/env python

import sys
from matchmaker.shmm_shmm_lib import *
from matchmaker.alignments_lib import *

def compute_q_combined(pair, alignment_file):
  test_aln = Alignment()
  test_aln.readFromFile(alignment_file)
  ref_aln = Alignment()
  ref_aln.readFromFile(get_reference_alignment_path(pair))
  return Qcombined(test_aln, ref_aln)
  
def main():
  if len(sys.argv) < 3:
    print "Usage: %s pair alignment_file" % sys.argv[0]
    sys.exit(0)

  pair = sys.argv[1]
  alignment_file = sys.argv[2]
  print compute_q_combined(pair, alignment_file)

if __name__ == '__main__':
  main()
