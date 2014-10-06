#! /usr/bin/env python
# Get the YL scores for all candidate matchmaker alignments for a pair 

import os, sys, glob
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

outfilepath = os.path.join(align_dir(seedX, seedY), "%s.yl" % (pair))
outfile = open(outfilepath, "w")

alist = matchmaker_seed_alignment_filenames(seedX, seedY,
    use_nr=(alignment_source=="nr"),
    use_consensus=(alignment_source=="consensus"),
    use_best=(alignment_source=="best"))

for afile in alist:
  afile_basename = os.path.basename(afile)
  (seed1, hmm1, seed2, hmm2, seeds_or_ambassadors) = \
    split_alignment_file_basename(afile_basename)
  if seeds_or_ambassadors == "YL":
    if seed1 != seedX:
      raise ShmmShmmLibError, \
       "seeds in filename do not match those on command line -- " + \
       "seedX = %s, seed1=%s (should match)" % (seedX, seed1)
    else:
      hmmX = hmm1
      hmmY = hmm2

    yl_path = yl_filename(seedX, seedY, hmmX, hmmY)
    if yl_path and os.path.exists(yl_path):
      yl_score = read_yl_score(yl_path)
      if yl_score:
	print >> outfile,  "%s %.5f" % (afile_basename, yl_score)

outfile.close()
