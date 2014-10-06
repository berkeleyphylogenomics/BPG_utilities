#!/usr/bin/python

"""
Written by Ruchira Datta, June 2008. Sjolander Lab.

Input: two seeds.
  Assumes that:
    (a) alignments have been generated
        following the matchmaker alignment protocol
        and stored in the protocol's alignment dir tree 
    (b) secondary structure predictons or assessments have been
        generated and stored in the protocol's "single_seqs"
        directory tree

What script does:
  Derives an alignment of secondary structure strings from
  each pairwise sequence alignment, and assesses the quality
  of that alignment.

Output: CSV file showing quality scores for each alignment,
    one alignment per record
"""

import os, sys, fasta, string
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 3:
  print "Usage: %s <seedX> <seedY>" % sys.argv[0]
  sys.exit(0)

seedX = sys.argv[1]
seedY = sys.argv[2]

def read_ss(seed):
  """ Read sec struc of seed from file in defined location.

      Return it in a string.
  """
  ss_file = os.path.join(single_seqs_dir(), "%s/%s.ss2" % (seed, seed))
  ss = {}
  f = open(ss_file)
  lines = f.readlines()
  f.close()
  for line in lines:
    line = line.rstrip()
    if len(line) > 0 and line[0] != '#':
      fields = line.split()
      ss[string.atoi(fields[0])] = fields[2]
  return ss

# read the secondary structure for each seed
X_ss = read_ss(seedX)
Y_ss = read_ss(seedY)

basedir = align_dir(seedX, seedY)
  
# initialize output CSV file: open and write header
out_file = "%s_%s_ss2.csv" % (seedX, seedY)
outfp = open(out_file, "w")
outfp.write("Alignment,NumAlignedPairs,")
outfp.write("NumAlignedPairsWithSameSecondaryStructure,")
outfp.write("FractionOfAlignedPairsHavingSameSecondaryStructure\n")

# this function is run on every directory containing alignments
def check_ss(arg, dirname, names):
  # break off the last three components of the dir pathname into
  # align, hmmX, and hmmY
  fulldirname = os.path.join(basedir, dirname)
  rest, hmmY = os.path.split(dirname)
  rest, hmmX = os.path.split(rest)
  rest, align = os.path.split(rest)
  print "===="
  print hmmY
  print hmmX
  print align
  idx = align.find('_align')

  # we want only SHMM-SHMM alignments (why?); therefore we want alignments
  # that are two levels down from the <seedX>_<seedY>_align dir
  # (i.e. where the third from last pathname component contains _align).
  if idx >= 0:
    # construct the alignment filename
    truename = "%s_%s_%s_%s_YL_original.afa" % (seedX, seedY, hmmX, hmmY)
    print truename
    # if this file exists in this dir ...
    if truename in names:
      print "FOUND!"
      # build the full pathname of the alignment file
      fulltruename = os.path.join(fulldirname, truename)
      # read the sequences from the file
      # into a 4-element list: header, seq, header, seq
      aligned_seqs = fasta.ReadSequencesList(fulltruename)
      aligned_X_seq = ""
      aligned_Y_seq = ""
      # figure out which seq is for seedX and which for seedY
      for header, aligned_seq in aligned_seqs:
        if len(header) >= len(seedX) and header[0:len(seedX)] == seedX:
          aligned_X_seq = aligned_seq
        elif len(header) >= len(seedY) and header[0:len(seedY)] == seedY:
          aligned_Y_seq = aligned_seq
        else: 
          print "Error: found extra sequence header %s in %s" % (header, 
                                                                  fulltruename)
      # if the lengths of the aligned sequences are identical
      # (which they should be), assess the quality of the
      # corresponding alignment of their secondary structure strings
      # and write this info in a single record of the output CSV file
      if len(aligned_X_seq) == len(aligned_Y_seq):
        i = 0
        j = 0
        num_aligned_pairs = 0
        num_ss_aligned_pairs = 0
        for k in range(len(aligned_X_seq)):
          cX = aligned_X_seq[k]
          cY = aligned_Y_seq[k]
          if cX.isupper() and cY.isupper():
            num_aligned_pairs = num_aligned_pairs + 1
            if X_ss[i + 1] == Y_ss[j + 1]:
              num_ss_aligned_pairs = num_ss_aligned_pairs + 1
          if cX.isalpha():
            i = i + 1
          if cY.isalpha():
            j = j + 1
        outfp.write("%s," % truename)
        if num_aligned_pairs == 0:
          outfp.write("0,0,0\n")
        else:
          frac = num_ss_aligned_pairs / float(num_aligned_pairs)
          outfp.write("%d,%d,%3f\n" % (num_aligned_pairs, num_ss_aligned_pairs,
                                        frac))
      else:
        print "Error: lengths of aligned sequences in %s differ" % fulltruename
        print "%s, length %s" % (seedX, len(aligned_X_seq))
        print aligned_X_seq
        print "%s, length %s" % (seedY, len(aligned_Y_seq))
        print aligned_Y_seq

os.path.walk(basedir, check_ss, None)
outfp.close()
