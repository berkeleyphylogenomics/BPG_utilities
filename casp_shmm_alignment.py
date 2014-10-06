#! /usr/bin/env python
# coding: iso-8859-15

""" casp_shmm_alignment.py       Terry Farrah    March 2008

Input: a target and template sequence, and a collection
   of ghmms/hmms for each
Output: a small set of proposed sequence alignments
  
Method:
-- do all vs. all alignments of HMMs and extract
     an alignment of the original pair from each HMM alignment
-- for each domain:
     * reduce alignment set via clustering and consensus
     * create a homology model based on each alignment
        and get DOPE score
     * output few alignments with highest DOPE scores
        and also output scores of alignments against reference alignment.

"""

import sys, os, glob, pdb, socket
import time
import Numeric
from matchmaker.shmm_shmm_lib import *
from matchmaker.build_model import *
from optparse import OptionParser

class TestShmmAlignment(Exception):
   pass

#====================
# Process command line
#====================
parser = OptionParser()
parser.add_option("-p","--pair",dest="pair",metavar="PAIR",
  help="Input seed pair (for example, 1bqsA_1mnc)")
parser.add_option("-a","--align",
  dest="do_align",
  action="store_true",
  default=False,
  help="Generate alignments from HMMs")
parser.add_option("-n","--nr",
  dest="do_nr",
  action="store_true",
  default=False,
  help="Make nonredundant list of alignments")
parser.add_option("-r","--rank",
  dest="do_rank",
  action="store_true",
  default=False,
  help="Rank alignments")
parser.add_option("-c","--cluster",
  dest="do_cluster",
  action="store_true",
  default=False,
  help="Cluster alignments")
parser.add_option("-o","--nonmatchmaker",
  dest="do_nonmatchmaker",
  action="store_true",
  default=False,
  help="Create alignments using non-matchmaker methods")
parser.add_option("-m","--model",
  dest="max_models",
  metavar="N",
  type="int",
  default=0,
  help="Create models for alignments, using up to the top " + \
       "N matchmaker alignments plus a few 'special' nonmatchmaker. " + \
       "When using this option, cannot submit " + \
       "this script to queue.")
(options, args) = parser.parse_args()

# check that necessary options are given, and that values are valid
# assign option values to variables
if not options.pair:
   parser.error("Seed pair required.")
pair = options.pair
do_align = options.do_align
do_nr = options.do_nr
do_rank = options.do_rank
do_cluster = options.do_cluster
do_nonmatchmaker = options.do_nonmatchmaker
max_models = options.max_models
seed_pair = pair.split("_")
seedX = seed_pair[0]
seedY = seed_pair[1]

# echo option choices/defaults
print "Pair = %s" %  pair

#======================
# Do the work
#======================

# NOTE:
# if we are using a new set of parameters, we need to
# create new dir like top500_unique90_cropped_gathered_homologs
# and create symbolic link "current_method" to it.

os.chdir(shmm_root_dir())

def do_step(cmd):
   cwd = os.getcwd()
   print >> sys.stderr, "cwd =", cwd
   print >> sys.stderr, "Running:", cmd
   os.system(cmd)

# required before running pipeline:
# -- pdb file for template
# -- seq file for target and template


#--------------------------------------------------
#  Derive pairwise alignments
#--------------------------------------------------

# Align seed seq, GHMM, GHMM consensus, and all SHMMs and
# SHMM consensuses for each seed with
# those for other seed in all pairwise combos
# and derive a seed-seed pairwise alignment from each
# input: *.mod from above, seed seqs, and ref seqs

if do_align:
   align_cmd = "align_originals_via_shmms_and_score.py"
#   cmd = "%s -k --use_kerfw05_hmms -p 15 --pairs %s" % (align_cmd, pair)
#   do_step(cmd)
   cmd = "%s -k --use_kerfw05_hmms -p 20 --pairs %s " % (align_cmd, pair) \
        + "--profile_profile_scoring_function hhalignsamss"
   do_step(cmd)
#   cmd = "%s -k --use_kerfinfoshare_hmms -p 15 --pairs %s" % (align_cmd, pair)
#   do_step(cmd)
#   cmd = "%s -k --use_kerfinfoshare_hmms -p 20 --pairs %s" % (align_cmd, pair)
#   do_step(cmd)
   #cmd = "%s --use_sciphy_subfams --pairs %s" % (align_cmd, pair)
   #do_step(cmd)

   # Count all alignments created above
   cmd = "count_pairwise_alignments_for_pair.py %s" % pair
   do_step(cmd)

if do_nr:
   # Uniquify the alignments
   cmd = "print_nr_alignments_with_dict.py %s %s" % (seedX, seedY)
   do_step(cmd)


#--------------------------------------------------
#  Cluster alignments
#--------------------------------------------------

# Reduce size of alignment set by clustering/consensus
# Runs in 1-18 minutes on bpg compute nodes
if do_cluster:
   cmd = "reduce_alignment_set.py %s" % pair
   do_step(cmd)

#--------------------------------------------------
#  Create non-matchmaker alignments
#--------------------------------------------------

if do_nonmatchmaker:
  cmd = "casp_non_matchmaker_align.py -p %s" % pair
  do_step(cmd)

#--------------------------------------------------
#  Rank all alignments according to secondary structure agreement & YL score
#--------------------------------------------------

# Perl scripts below contributed by Sriram June 2008
# Expects the following:
#  -- files named YL.out throughout <align_dir_for_pair> hierarchy
#   (these are created dring the alignment process above)
#  -- .ss2 file for each protein in single_seq_dir()/<pair>
#   (this is created in the single-seed phase of the pipeline)

# 6/25/08: if clustering was not done or not successful,
#   list of nonredundant alignments will be ranked.
#  If uniquifying was not done, all alignments will be ranked.
#  This is implemented with the "-a best" option in 2 commands
#   below.

if do_rank:

   align_dir_for_pair = align_dir(seedX, seedY)
   os.chdir(align_dir_for_pair)
   # Remove any old rank files sitting around.
   # If we do this, it doesn't make sense to check for existance
   # of rank files before doing computations below.
   # Ranking is fast, and if user chooses rank option, she
   # probably wants the ranking to be done even if done before.
   cmd = "remove_rank_files.py %s %s" % (seedX, seedY)
   do_step(cmd)

   # Get the Yona-Levitt alignment quality score for each consensus alignment
   cmd = "topdir='%s'" % (align_dir_for_pair)
   do_step(cmd)
   if not os.path.exists("%s.yl" % (pair)):
     # gets all files *YL.out in directory tree align_dir_for_pair
     cmd = "get_yl_per_pair.py -p %s -a best" % (pair)
     do_step(cmd)

   # Compute the secondary structure similarity score for
   # each alignment
   if not os.path.exists("%s_ss2.csv" % (pair)):
     cmd = "compute_ss_score.py -p %s -a best" %  (pair)
     do_step(cmd)
   if not os.path.exists("%s.ss2" % (pair)):
     cmd = "cleanss2.py %s_ss2.csv > %s.ss2" % (pair, pair)
     do_step(cmd)

   # Apply Sriram's ranking method to this data
   if not os.path.exists("%s.features" % (pair)):
     cmd = "integrate_features.pl %s %s > %s.features" % (seedX, seedY, pair)
     do_step(cmd)
   if not os.path.exists("%s.features.ext" % (pair)):
     cmd = "extract_relevant_features.pl %s.features > %s.features.ext" \
       % (pair, pair)
     do_step(cmd)
   if not os.path.exists("%s.rank" % (pair)):
     cmd = "rank.pl %s.features.ext | sort -n -r -t $'\t' -k 2 > %s.rank" % \
        (pair, pair)
     do_step(cmd)

   # Ranking is now in <align_dir_for_pair>/<pair>.rank


#--------------------------------------------------
#  Select best alignments using modeling & scoring
#--------------------------------------------------

if max_models:
   
   # all of the steps below work on the top max_models ranked alignments
   # plus 5 additional alignments

   # Mask alignments if necessary.
   cmd = "mask_alignments.py -p %s -n %d" % (pair, max_models)
   do_step(cmd)

   # Copy PDB file to central repository
   # Global variables pdb_dir, pdb_prefix, and pdb_ext are defined
   # in matchmaker.build_model
   pdb_source_filename = seedX[0:4] + ".pdb"
   pdb_source_path = os.path.join(single_seqs_dir(), seedX,
        pdb_source_filename)
   pdb_destination_filename = pdb_prefix + seedX[0:4] + "." + pdb_ext
   pdb_destination_path = os.path.join(pdb_dir, pdb_destination_filename)
   cmd = "cp %s %s" % (pdb_source_path, pdb_destination_path)
   do_step(cmd)

   # Create rough model for each algnmt w/coarse refinement
   #  & obtain TSVmod scores
   cmd = "build_and_score_models_for_pair.py %s %d" % (pair, max_models)
   do_step(cmd)
