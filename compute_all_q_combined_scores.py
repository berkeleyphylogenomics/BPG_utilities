#!/usr/bin/env python

import os, sys
from matchmaker.shmm_shmm_lib import *
from matchmaker.compute_q_combined_score import compute_q_combined

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <pair>" % sys.argv[0]
    sys.exit(0)

  pair = sys.argv[1]
  seedX, seedY = pair.split('_')
  for alignment_file in matchmaker_seed_alignment_filenames(seedX, seedY, 
                                                             use_nr=True):
    if os.path.exists(alignment_file):
      qcombined = compute_q_combined(pair, alignment_file)
      print "%s,%g" % (os.path.split(alignment_file)[1], qcombined)

if __name__ == '__main__':
  main()
