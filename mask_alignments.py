#!/bin/env python

"""
Looks for a mask file for template sequence, indicating that certain
residues in that sequence do not have coordinates in the corresponding
PDB file. Applies that mask to all alignments for target/template pair.
"""
import os,sys,glob
import BPG_common.fasta
from matchmaker.shmm_shmm_lib import *
from optparse import OptionParser

max_ranked = 200

def mask_alignment_file(alignment_filename, mask, pair):
  if alignment_filename.startswith("/"):
    # we have full pathname for nonmatchmaker alignment
    fulldirname = os.path.dirname(alignment_filename)
    alignment_filename = os.path.basename(alignment_filename)

  else:
    # otherwise, this is a matchmaker alignment filename
    fulldirname = alignment_dir(alignment_filename)

  os.chdir(fulldirname)
  unmasked_alignment_filename = os.path.splitext(alignment_filename)[0] + \
      ".unmasked.a2m"
  # this file may already exist if we have masked before.
  if not os.path.exists(unmasked_alignment_filename):
    # rename .afa file
    cmd = "mv %s %s" % (alignment_filename, unmasked_alignment_filename)
    os.system(cmd)
  print unmasked_alignment_filename
  (seedX, aligned_X_seq), (seedY, aligned_Y_seq) = \
    read_alignment_file(unmasked_alignment_filename)
  write_alignment_file_compressed(mask_aligned_seq(aligned_X_seq, mask),
       aligned_Y_seq, alignment_filename)

def main():
  parser = OptionParser()

  parser.add_option("-a", "--alignment_source", dest="alignment_source",
	metavar="SOURCE",
	type="choice",
	choices=["ranked_spec", "ranked", "consensus", "nr", "all", "best"],
	default="ranked_spec",
  help="Specify which alignment source to use. " + \
	"Choices are: ranked_spec (top ranked + special, default), ranked (top ranked consensus alignments), consensus (consensus alignments after clustering), nr (non-redundant), all, best (smallest set available from consensus/nr/all)."
        )
  parser.add_option("-n", "--num_ranked", dest="num_ranked",
	metavar="N",
	type="int",
        default=200,
  help="maximum number of ranked alignments to use (default 200)")
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
  max_ranked = options.num_ranked
  alignment_source = options.alignment_source

  (template, target) = pair.split("_")

  mask = get_mask(template)
  # if there is no mask file, then no masking is needed.
  if mask == None: sys.exit(0)

  initial_dir = os.getcwd()

  # Get file listing consensus alignments in ranked order.
  # Want to generalize this so that we can use other kinds of input
  # (all alignments, uniquified alignments, unranked consensus alignments) --
  # perhaps via command line options
  if alignment_source == "ranked_spec":
    alignment_list = \
       get_top_ranked_and_special_alignment_filenames(pair, max_ranked)
  elif alignment_source == "ranked":
    alignment_list = get_top_ranked_alignment_filenames(pair, max_ranked)
  else:
    alignment_list = all_seed_alignment_filenames(seedX, seedY,
      use_nr=(alignment_source=="nr"),
      use_consensus=(alignment_source=="consensus"),
      use_best=(alignment_source=="best"))

  for alignment_filename in alignment_list:
    mask_alignment_file(alignment_filename, mask, pair)

  os.chdir(initial_dir)

if __name__ == "__main__":
  main()
