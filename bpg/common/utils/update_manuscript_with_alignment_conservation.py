#!/usr/bin/env python

import os, sys
from bpg.makebook.buildFamily import computeAlignmentConservation
from Bio import SeqIO

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <alignment_path>" %  sys.argv[0]
    sys.exit(0)

  alignment_path = sys.argv[1]
  basename = os.path.splitext(alignment_path)[0]

  handle = open(alignment_path, 'r')
  alignmentrecords = list(SeqIO.parse(handle, 'fasta'))
  handle.close()
  computeAlignmentConservation(basename, alignmentrecords)

if __name__ == '__main__':
  main()
