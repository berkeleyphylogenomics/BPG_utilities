#! /usr/bin/env python
# coding: iso-8859-15

"""
casp_create_hmms.py       Terry Farrah    May 2008

Input: a protein sequence -- if CASP target, must be in the CASP dir
  -- if template, must be in PDB
  (a sequence universe database must also be available)
Output: ghmms and shmms for the protein's family, in preparation for
  SHMM-SHMM alignment protocol
  
Method:
-- collect other family members for each seq, subdivide
   each family, and create an HMM for each subfamily

"""

import sys, os, glob, pdb, socket, stat
import time, cPickle
import Numeric
from matchmaker.shmm_shmm_lib import *
from optparse import OptionParser
from Bio import Seq, SeqIO

casp_target_dir = "/clusterfs/ohana/software/CASP8/targets"
SABmark_base_dir = "/clusterfs/ohana/external/SABmark/twi"
rcsb_website = "http://www.rcsb.org/pdb/files/"
pdb_db_name = "pdb_rcsb"

class CaspShmmAlignment(Exception):
   pass

#====================
# Process command line
#====================
parser = OptionParser()
parser.add_option("-S","--SABmark", dest="use_SABmark", action="store_true",
                    default=False,
                    help="Use sequence from SABmark reference alignment extended to full PDB sequence")
parser.add_option("-s","--seq_id",
  dest="seed",
  metavar="SEQ_ID",
  help="Input <seq_id>. If template, must be PDB identifier. If target, must be CASP8 target identifier, and seq must be in a .fa file in the CASP8 directory.")
parser.add_option("-t", "--type",
  dest="seed_type",
  metavar="SEED_TYPE",
  type="choice",
  choices=["target","template"],
  help = "choose target (structure unknown) or template (structure known)")
 
(options, args) = parser.parse_args()

# check that necessary options are given, and that values are valid
# assign option values to variables
if not options.seed:
   parser.error("Seed ID not provided.")
seed = options.seed
if not options.seed_type:
  parser.error("Option -t required.")
seed_type = options.seed_type

# echo option choices/defaults
print "Seq ID = %s" %  seed
print "Seed type = %s" % seed_type

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
   print  "cwd =", cwd
   print  "%s: running '%s'" % \
          (time.asctime(time.localtime()), cmd)
   os.system(cmd)

#----------------------------------------------------
#  Get sequence
#----------------------------------------------------

# Create directory for files related to the sequence

# Adjust seed names to be of this form:
# 4-char ID, followed by single-char chainID for template
if options.use_SABmark:
  f = open("/clusterfs/ohana/external/pdb_of_scop_1.65.pkl")
  pdb_of_scop = cPickle.load(f)
  f.close()
  if seed not in pdb_of_scop:
    print "%s not found in pdb_of_scop file" % seed
    sys.exit(1)
  original_seed = seed
  seed_without_chainID, chainID, start, end = pdb_of_scop[seed]
  seed_pdb_name = seed_without_chainID + "_" + chainID
  seed = seed_without_chainID + chainID
else:
  if seed_type == "template":
    if len(seed) == 4:
      seed_without_chainID = seed
      chainID = "A"
      seed_pdb_name = seed + "_" + chainID
      seed = seed + chainID
    elif len(seed) == 5:
      chainID = seed[4]
      seed_pdb_name = seed[0:4] + "_" + chainID
      seed_without_chainID = seed[0:4]
  if seed_type == "target":
    target_name = seed
    seed = 't' + target_name[2:5]

# Make a directory for data related to a single sequence
single_seq_dir = os.path.join(single_seqs_dir(), seed)
make_dir_exist(single_seq_dir)
os.chdir(single_seq_dir)

if options.use_SABmark:
  # If options.use_SABmark is enabled, get sequence from the SABmark full
  # alignment (i.e., the SABmark alignment extended to the full PDB sequence;
  # this will also include any residues from the SABmark alignment that are not
  # in the full PDB sequence, perhaps due to a break in the PDB chain).
  cmd = "grep %s %s/group*/group.summary" % (original_seed, SABmark_base_dir)
  status, output = commands.getstatusoutput(cmd)
  try:
    group = output.split('/')[6]
    alignment_file = glob.glob("%s/%s/reference/*%s*full.fasta" 
                              % (SABmark_base_dir, group, original_seed))[0]
    f = open(alignment_file, "rU")
    records = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    f.close()
    record = records[original_seed]
    unaligned_seq = record.seq.tostring().replace('-', '')
    seed_seq_filename = os.path.join(single_seq_dir, "%s.fa" % seed)
    f = open(seed_seq_filename, "w")
    f.write(">%s\n" % seed_pdb_name)
    f.write("%s\n" % unaligned_seq)
    f.close()
  except KeyError:
    print "Full SABmark reference alignment for %s not found" % seed
    sys.exit(1)

else:
  # Get sequence and put it in that dir, changing header
  # so that it just contains seeed name
  # If template, get sequence from PDB via internet.
  # If target, get sequence from CASP data directory
  # (currently, needs to be put there manually)
  if seed_type == "template":
    seed_seq_filename = "%s.fa" % (seed)

    """ 6/13/08: we are doing the below manually now
        because we may need to add/subtract residues
        from the sequence file, and because we can't
        download web files from the nodes
    if not os.path.exists(seed_seq_filename) or \
         os.path.getsize(seed_seq_filename) == 0:
      # change header line so that it just contains the seed name
      cmd = "fastacmd -d %s -s %s | sed '/>/c\>%s'> %s"  \
       % (pdb_db_name, seed_pdb_name, seed, seed_seq_filename)
      do_step(cmd)


    # Get PDB file from internet
    # 6/13/08: wget and curl do not work on nodes
    pdb_filename = "%s.pdb" %  (seed_without_chainID)
    pdb_try = 0
    max_pdb_try = 5
    while not os.path.exists(pdb_filename) and pdb_try < max_pdb_try:
      cmd = "curl ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz > %s.gz" \
         % (seed[1:3], seed_without_chainID, pdb_filename)
      # wget does not exist on most ohana nodes
      #cmd = "/usr/bin/wget %s" %  (os.path.join(rcsb_website, pdb_filename))
      if pdb_try > 0: os.system("sleep 30")
      do_step(cmd)
      pdb_try = pdb_try + 1
    if not os.path.exists("%s.gz" % pdb_filename):
      print >> sys.stderr, "Could not download %s file after %d tries" \
        % (pdb_filename, max_pdb_try)
    cmd = "gunzip %s.gz" % (pdb_filename)
    do_step(cmd)

    """

    # Remove records for chains other than the one we want,
    # and store remaining PDB records in .ent file for dssp
    cmd = "get_pdb_single_chain.py %s < %s.pdb > %s.ent" % \
       (chainID, seed_without_chainID, seed)
    do_step(cmd)

  # if target
  else:
    seed_seq_filename = "%s.fa" % (seed)
    if not os.path.exists(seed_seq_filename):
      target_seq_filename = \
          os.path.join(casp_target_dir, target_name, "%s.fa" % target_name)
      if not os.path.exists(target_seq_filename):
        print >> sys.stderr, "Target seq file, %s, does not exist." %  \
            (target_seq_filename)
        sys.exit(1)
      # change header line so that it just contains the seed name
      cmd = "cat %s | sed '/>/c\>%s'> %s"  % (target_seq_filename, seed,
                seed_seq_filename)
      do_step(cmd)



#----------------------------------------------------
#  Gather family for each seed; create and refine MSA
#----------------------------------------------------

# collect homologs for seed using FlowerPower
# input: seed sequence in <shmm_root>/single_seqs/<seed>/<seed>.fa
# output: family MSA in <shmm_root>/flowerpower_seedmaster/<seed>/final.a2m
fp_dir = os.path.join(shmm_root_dir(), "flowerpower_seedmaster", seed)
make_dir_exist(fp_dir)
outfile = os.path.join(fp_dir, "final.a2m")
if not os.path.exists(outfile) or os.stat(outfile)[stat.ST_SIZE] == 0:
  os.chdir(fp_dir)
  infile = "%s/%s/%s.fa" % (single_seqs_dir(), seed, seed)
  if not os.path.exists(infile):
    print >> sys.stderr, "Seed sequence file %s does not exist." % (infile)
    sys.exit(0)
  cmd = "build_family %s" % (infile)
  do_step(cmd)
else:
  print "FlowerPower already run on %s" % seed

# If FP ran successfully, delete workspace directory
#if os.path.exists(outfile):
#  fp_work_dir = os.path.join(fp_dir, "fp-workspace")
#  cmd = "rm -r -f %s" % fp_work_dir
#  do_step(cmd)

os.chdir(shmm_root_dir())

# link the final FlowerPower MSA to another filename and location
# output: <shmm_root>/gathered_homologs/<seed>/<seed>_homologs.a2m
outfile = os.path.join(shmm_root_dir(), "gathered_homologs", seed,
    "%s_homologs.a2m" % (seed))
if not os.path.exists(outfile):
  cmd = "link_gathered_homologs.py %s" % seed
  do_step(cmd)

# Crop MSA to seed
# output:
# <shmm_root>/cropped_gathered_homologs/<seed>/<seed>_homologs_cropped_to_seed.a2m"
outfile = os.path.join(shmm_root_dir(), "cropped_gathered_homologs", seed,
   "%s_homologs_cropped_to_seed.a2m" % (seed))
if not os.path.exists(outfile):
  cmd = "crop_to_seed.py %s" % seed
  do_step(cmd)

# Make MSA nonredundant at 90% using uniqueseq
# output: # <shmm_root>/nr100_cropped_gathered_homologs/<seed>/<seed>_cropped_nr100.a2m
# 07/20/09: Stopped using unique
infile = outfile
outdir = os.path.join(shmm_root_dir(), "nr100_cropped_gathered_homologs", seed)
make_dir_exist(outdir)
outfile = os.path.join(outdir, "%s_homologs_cropped_nr100.a2m" % (seed))
if not os.path.exists(outfile):
  current_directory = os.getcwd()
  os.chdir(outdir)
  cmd = "make_nr_at_100_with_dict.py %s_homologs_cropped %s" \
        % (seed, infile)
  do_step(cmd)
  cmd = "ln -s %s_homologs_cropped_nr100.fa %s_homologs_cropped_nr100.a2m" \
          % (seed, seed)
  do_step(cmd)
  os.chdir(current_directory)

# Take top N seqs according to %ID with seed
# output: <shmm_root>/top<N>_unique90_cropped_gathered_homologs/<seed>/<seed>.a2m
#    = <cur_dir>/<seed>/<seed>.a2m  (final MSA for family)
# 4/20/08: this does not run on the nodes because it calls site-specific
#   python libs
# 4/19/08: we would like to leave this step out; in order to do so,
#   we need to change definition of cur_dir and rename MSAs from previous step
#   "<seed>.a2m"
# 5/13/08: step removed
# 6/9/08: step reinstituted, with 6000 instead of 500
#   (effectively takes all seqs; purpose of step is to add seed back in)
outdir = os.path.join(shmm_root_dir(),
    "top6000_nr100_cropped_gathered_homologs")
outfile = os.path.join(outdir, seed, "%s.a2m" % (seed))
if not os.path.exists(outfile):
  cmd = "top_by_pwid_with_seed.py 6000 %s" % seed
  do_step(cmd)

# Make symbolic link current_method point here;
# succeeding scripts will place results here and look for the
# final MSA here.
# 6/12/08: this caused trouble today. Somehow, current_method
#  got created as a real dir. Could not then be removed in
#  command below. Link was made to *inside* existing dir
#  current_method. Messed up placement of subsequent output files.
#cmd = "rm current_method; ln -s %s current_method" % (outdir)
#do_step(cmd)


#----------------------------------------------------
#  Divide each family into subfamilies
#----------------------------------------------------

# in comments below, <cur_dir> =
#   <shmm_root>/top<N>_unique90_cropped_gathered_homologs =
#   <shmm_root>/current_method   (via symbolic link)

# First, create GHMM for family
# input: <cur_dir>/<seed>/<seed>.a2m (MSA for family)
# output: <cur_dir>/<seed>/<seed>.mod
ghmmfile = os.path.join(cur_dir(), seed, "%s.mod" % seed)
if not os.path.exists(ghmmfile):
  cmd = "build_ghmm.py %s" % seed
  do_step(cmd)

# Convert ghmm to ascii
# input: <cur_dir>/<seed>/<seed>.mod
# output: <cur_dir>/<seed>/ascii.ghmm/<seed>.mod
#         <cur_dir>/<seed>/kerf<N>/info-share/<seed>.mod (sym link to ascii GHMM)
ascii_ghmmfile = os.path.join(cur_dir(), seed, "ascii_ghmm", "%s.mod" % seed)
if not os.path.exists(ascii_ghmmfile):
  cmd = "make_ascii_ghmm.py %s" % seed
  do_step(cmd)


#--- Method #1: phylip-kerf -----------

# Create neighbor-joining phylogenetic tree for MSA using Phylip
# input: <cur_dir>/<seed>/<seed>.a2m (MSA for family)
# output: <cur_dir>/<seed>/<seed>.nj
outfile = os.path.join(cur_dir(), seed, "%s.nj" % seed)
if not os.path.exists(outfile):
  cmd = "make_nj_tree.py %s" % seed
  do_step(cmd)

# Create subfamilies by cutting trees at 15% and 20% ID using kerf
# output: <cur_dir>/<seed>/kerf<N>/<seed>.kerf<N>.sf<m>.{fa,tre}
#                    (MSA and tree for each subfam)
#         <cur_dir>/<seed>/<seed>.kerf<N>.sf<m>.fa (consensus for subfams)
#                      where m goes from 0 to M-1, M = # subfams
#         <cur_dir>/<seed>/<seed>.ghmm.fa (?) (consensus for GHMM)
# Tests for whether kerf has been run already are done within do_kerf.py
#cmd = "do_kerf.py 15 %s" % seed
#do_step(cmd)
cmd = "do_kerf.py 20 %s" % seed
do_step(cmd)


#--- Method #2: sci-phy (creates SHMMs, too)  -----------

# Create phylogenenetic trees and generate subfamilies/SHMMs using SCI-PHY
# input: <cur_dir>/<seed>/<seed>.a2m   (MSA for family)
#        <cur_dir>/<seed>/ascii.ghmm/<seed>.mod  (GHMM for family)
# output: <cur_dir>/<seed>/sciphy/<seed><sp_serial_num>.mod  (SHMMs)
# Test for whether this step has already been run is executed within
#  make_sciphy_subfams.py
#cmd = "make_sciphy_subfams.py %s" % seed
#do_step(cmd)


#--------------------------------------------------
#  Create SHMMs for kerf subfamilies 
#--------------------------------------------------

# Method #1: create SHMMs from kerf subfams using w0.5
# input: <cur_dir>/<seed>/kerf<N>/<seed>.kerf<N>.sf<m>.fa
# output: <cur_dir>/<seed>/kerf<N>/w05/<seed>.w05kerf<N>.sf<m>.mod
#           where m goes from 0 to M-1
# Test for whether this step has already been run is executed within
#  build_w05_kerf_shmms.py
#cmd = "build_w05_kerf_shmms.py 15 %s" % seed
#do_step(cmd)
cmd = "build_w05_kerf_shmms.py 20 %s" % seed
do_step(cmd)


# Method #2: create SHMMs from kerf subfams using sci-phy / info-sharing
# input: <cur_dir>/<seed>/kerf<n>/kerf.seqs
#   where n = all kerf %IDs previously generated
# output: <cur_dir>/<seed>/kerf<n>/sciphy/<seed>.infokerf<n>.<m>.mod
# Test for whether this step has already been run is executed within
#  build_infoshare_kerf_shmms.py
#cmd = "build_infoshare_kerf_shmms.py 15 %s" % seed
#do_step(cmd)
#cmd = "build_infoshare_kerf_shmms.py 20 %s" % seed
#do_step(cmd)

#=============================================
#  Create secondary structure file in prep for ranking
#=============================================

# Step for target requires HMM, so this code segment is
# placed after HMM creation.
os.chdir(single_seq_dir)
if seed_type == "template":
  cmd = "dsspcmbi %s.ent > %s.dssp" % (seed, seed)
  do_step(cmd)
  seq_file = "%s.fa" % seed
  mask_file = "%s.mask" % seed
  if os.path.exists(mask_file): dssp_seq_file = mask_file
  else: dssp_seq_file = seq_file
  cmd = "dssp_to_threestate.py -i %s.dssp -s %s -o %s.ss2" % (seed,
        dssp_seq_file, seed)
  do_step(cmd)
elif seed_type == "target":
  # may want to use shmm instead.
  cmd = "run_psipred_with_hmm.py %s %s" % \
    (seed, ghmmfile)
  do_step(cmd)

#=============================================
#  Report which GHMMs and SHMMs were created
#=============================================

print "================================================="
print "%s complete for %s" % (sys.argv[0], seed)
print "================================================="
print "The following subfamilies/SHMMs exist for %s:" % (seed)
cmd = "report_hmms_for_seed.py %s" % (seed)
do_step(cmd)
