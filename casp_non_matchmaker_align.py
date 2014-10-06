#! /usr/bin/env python
# coding: iso-8859-15

""" casp_nonshmm_alignment.py       Terry Farrah    June 2008

Input: a target and template sequence, and homologs, MSAs, and
   ghmms for each, stored in the standard locations for matchmaker
Output: sequence alignments using various non-matchmaker methods
  
Methods:
  seed-prof (generated during matchmaker)
    -- seq-prof
    -- prof-seq
  prof-prof
    -- HHalign
    -- compass
    -- profsim (generated during matchmaker)
  multipw
    -- amap-multipw
    -- clustalw-multipw
    -- mafft-multipw
    -- muscle-multipw
    -- proda-multipw

"""

import sys, os, glob, pdb, socket
import time
import Numeric
from matchmaker.shmm_shmm_lib import *
from optparse import OptionParser

class NonShmmAlignment(Exception):
   pass

#====================
# Process command line
#====================
parser = OptionParser()
parser.add_option("-p","--pair",dest="pair",metavar="PAIR",
  help="Input seed pair (for example, 1bqsA_1mnc)")

(options, args) = parser.parse_args()

# check that necessary options are given, and that values are valid
# assign option values to variables
if not options.pair:
   parser.error("Seed pair required.")
pair = options.pair
seed_pair = pair.split("_")
seedX = seed_pair[0]
seedY = seed_pair[1]

# echo option choices/defaults
print "Pair = %s" %  pair

#======================
# Do the work
#======================

def do_step(cmd):
   cwd = os.getcwd()
   print >> sys.stderr, "cwd =", cwd
   print >> sys.stderr, "Running:", cmd
   os.system(cmd)


#--------------------------------------------------
#  HHalign
#--------------------------------------------------

seqX_pathname = os.path.join(single_seqs_dir(), seedX, "%s.fa" % seedX)
seqX_filename = os.path.basename(seqX_pathname)
seqY_pathname = os.path.join(single_seqs_dir(), seedY, "%s.fa" % seedY)
seqY_filename = os.path.basename(seqY_pathname)
ghmmX_pathname = ghmm_file_of_seed(seedX)
ghmmX_filename = os.path.basename(ghmmX_pathname)
ghmmY_pathname = ghmm_file_of_seed(seedY)
ghmmY_filename = os.path.basename(ghmmY_pathname)

# Create directory for HHalign and copy data into it,
# because SAM is picky about this.
hhalign_dir = os.path.join(shmm_root_dir(), "hhalign", pair)
make_dir_exist(hhalign_dir)
os.chdir(hhalign_dir)
cmd = "cp %s ." % ghmmX_pathname
os.system(cmd)
cmd = "cp %s ." % ghmmY_pathname
os.system(cmd)
cmd = "cp %s ." % seqX_pathname
os.system(cmd)
cmd = "cp %s ." % seqY_pathname
os.system(cmd)

# Run HHalign without secondary structure prediction
#cmd = "hhalign_with_hmms.py --seq1=%s --seq2=%s --hmm1=%s --hmm2=%s -o %s_hhalign.afa >| hhalign.out" % \
#   (seqX_filename, seqY_filename, ghmmX_filename, ghmmY_filename, pair)
#do_step(cmd)
# Run HHalign with secondary structure prediction
cmd = "hhalign_with_hmms.py -s --seq1=%s --seq2=%s --hmm1=%s --hmm2=%s -o %s_hhalignss.afa >| hhalign.out" % \
   (seqX_filename, seqY_filename, ghmmX_filename, ghmmY_filename, pair)
do_step(cmd)
compress_alignment_file("%s_hhalignss.afa" % pair)

#cmd = "cat %s.hhalign" % pair

#--------------------------------------------------
#  multi-pw
#--------------------------------------------------

# prepare for multi-pw
cmd = "concat_for_multipw.py %s %s" % (seedX, seedY)
do_step(cmd)

# amap needs to be installed on ohana
#cmd = "amap-multipw.py %s %s" % (seedX,  seedY)
#do_step(cmd)

cmd = "clustalw-multipw.py %s %s" % (seedX,  seedY)
do_step(cmd)

# mafft needs to be installed on ohana
#cmd = "mafft-multipw.py %s %s" % (seedX,  seedY)
#do_step(cmd)

cmd = "muscle-multipw.py %s %s" % (seedX,  seedY)
do_step(cmd)

# proda needs to be installed on ohana
#cmd = "proda-multipw.py %s %s" % (seedX,  seedY)
#do_step(cmd)
