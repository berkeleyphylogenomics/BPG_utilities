#! /usr/bin/env python

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

import os, sys, string
from BPG_common.fasta import *
import glob
from matchmaker.shmm_shmm_lib import *
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-a", "--alignment_source", dest="alignment_source",
      metavar="SOURCE",
      type="choice",
      choices=["consensus", "nr", "all", "best"],
      default="consensus",
help="Specify which alignment source to use: consensus (from clustering), non-redundant, all, or best available (least redundant)")
parser.add_option("-p", "--pair", dest="pair",
      metavar="SEEDPAIR",
      type="string",
help="seed pair to use, of form template_target")

(options, args) = parser.parse_args()

# check that necessary options are given, and that values are valid
# assign option values to variables
if not options.pair:
   parser.error("Option -p required")
pair = options.pair
(seedX,  seedY) = pair.split("_")
alignment_source = options.alignment_source

initial_dir = os.getcwd()
os.chdir(align_dir(seedX,seedY))

# read the secondary structure for each seed into a dictionary
X_ss = read_ss(seedX)
Y_ss = read_ss(seedY)

# initialize output CSV file: open and write header
out_file = "%s_%s_ss2.csv" % (seedX, seedY)
outfp = open(out_file, "w")
outfp.write("Alignment,NumAlignedPairs,")
outfp.write("NumAlignedPairsWithSameSecondaryStructure,")
outfp.write("FractionOfAlignedPairsHavingSameSecondaryStructure\n")

def check_ss(alignment_filename, pair, outfp, matchmaker):
      # read the aligned sequences from the file
      # into a 4-element list: header, seq, header, seq
      """
      # 7/15/08:  replaced with function to read possibly compressed file
      aligned_seqs = ReadSequencesList(alignment_filename)
      """
      aligned_seqs = read_alignment_file(alignment_filename)
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
                                                         alignment_filename)
      # if the lengths of the aligned sequences are identical
      # (which they should be), compute the percent ID of the
      # corresponding alignment of their secondary structure strings
      # and write this info in a single record of the output CSV file

      if len(aligned_X_seq) == len(aligned_Y_seq):
        (num_aligned_pairs, num_ss_aligned_pairs, frac) = \
         sec_structure_identity(aligned_X_seq, aligned_Y_seq, X_ss, Y_ss)
          
        # if this is a matchmaker alignment, we only need to write
        # the basename
        if matchmaker:
          outfp.write("%s," % os.path.basename(alignment_filename))
        else:
          outfp.write("%s," % alignment_filename)
        if num_aligned_pairs == 0:
          outfp.write("0,0,0\n")
        else:
          outfp.write("%d,%d,%3f\n" % (num_aligned_pairs, num_ss_aligned_pairs,
                                        frac))
      else:
        print "Error: lengths of aligned sequences in %s differ" % \
                 alignment_filename
        print "%s, length %s" % (seedX, len(aligned_X_seq))
        print aligned_X_seq
        print "%s, length %s" % (seedY, len(aligned_Y_seq))
        print aligned_Y_seq

mm_alist = matchmaker_seed_alignment_filenames(seedX, seedY,
    use_nr=(alignment_source=="nr"),
    use_consensus=(alignment_source=="consensus"),
    use_best=(alignment_source=="best"))
for alignment_filename in mm_alist:
  check_ss(alignment_filename, pair, outfp, matchmaker=True)

non_mm_alist = non_matchmaker_seed_alignment_filenames(seedX, seedY,
    use_nr=(alignment_source=="nr"))
for alignment_filename in non_mm_alist:
  check_ss(alignment_filename, pair, outfp, matchmaker=False)

os.chdir(initial_dir)
