#!/usr/bin/python

import os, sys
from Bio import Seq, SeqIO
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 3:
  print "Usage: %s <seedX_id> <seedY_id>" % sys.argv[0]
  sys.exit(0)

seedX_id = sys.argv[1]
seedY_id = sys.argv[2]
pair_id = "%s_%s" % (seedX_id, seedY_id)

# FlowerPower has gathered homologs for each seed and created an alignment.
# These have been concatenated.
out_dir = "%s/proda-multipw/%s" % (results_dir(), pair_id)
make_dir_exist(out_dir)
os.chdir(out_dir)
os.system("ln -s ../../../concat_for_multipw/%s/%s_concatenated.a2m " \
  % (pair_id, pair_id) +
  "%s_concatenated.mfa" % pair_id)
# Align this concatenation of the two MSAs with ProDA.
out_name = "%s_prodamulti.afa" % pair_id
# 6/17/08: proda  needs to be installed on ohana
os.system("proda -fasta %s_concatenated.mfa > proda_out" % pair_id)
os.system("mv %s_concatenated.fasta %s" % (pair_id, out_name))
proda_out_path = os.path.join(out_dir, out_name)
# Extract the pairwise alignment of the seeds from the new MSA.
seedX_seq = ""
seedY_seq = ""
for record in SeqIO.parse(open(proda_out_path, "rU"), "fasta"):
  if record.id == seedX_id or record.id == seedX_id[0:-1]:
    seedX_seq = record.seq.tostring()
  elif record.id == seedY_id or record.id == seedY_id[0:-1]:
    seedY_seq = record.seq.tostring()
if seedX_seq == "":
  print "Couldn't find %s in the MSA, exiting\n" % seedX_id
  sys.exit(0)
if seedY_seq == "":
  print "Couldn't find %s in the MSA, exiting\n" % seedY_id
  sys.exit(0)
# Output the pairwise alignment
pw_alignment_path = os.path.join(out_dir, "%s_prodamultipw.afa" % pair_id)
fp = open(pw_alignment_path, "w")
fp.write(">%s\n" % seedX_id)
fp.write("%s\n" % seedX_seq)
fp.write(">%s\n" % seedY_id)
fp.write("%s\n" % seedY_seq)
fp.close()

# Directory where structural alignments are kept.

structural_alignment_dir = "/home/ruchira/SHMM-SHMM/pairs"

# If there is a structural alignment, generate CS score

reference_alignment_path = os.path.join(structural_alignment_dir, 
                                        pair_id, "%s.afa" % pair_id)
if os.path.exists(reference_alignment_path):
  cs_output_file = os.path.join(out_dir, "%s_lobster_CS.out" % pair_id)
  os.system("/home/ruchira/bin/lobster -compareMSA %s -ref %s >& %s " \
	    % ( pw_alignment_path, reference_alignment_path, 
		cs_output_file ))
  (sp_score, cs_score) = read_sp_cs_score(cs_output_file)
  output_file = os.path.join(out_dir, "%s_proda-multipw.csv" % pair_id)
  fp = open(output_file, "w")
  fp.write("pair,SP,CS\n")
  fp.write("%s,%0.3f,%0.3f\n" % (pair_id, sp_score, cs_score))
  fp.close()

