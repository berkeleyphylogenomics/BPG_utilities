#!/usr/bin/python

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

"""
June 17, 2008  Terry Farrah  Sjolander Lab

Part of the matchmaker (SHMM-SHMM) alignment pipeline.
Counts the pairwise matchmaker alignments between seed sequences
found in the alignment directory tree for the pair.

"""

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <pair>" % sys.argv[0]
    sys.exit(0)

  pair = sys.argv[1]
  (seedX, seedY) = pair.split("_")
  n = len(matchmaker_seed_alignment_filenames(seedX, seedY))
  print "%d pairwise matchmaker alignments found for %s" % (n, pair)



if __name__ == "__main__":
  main()
