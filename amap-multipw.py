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
amap_in_path = "%s/concat_for_multipw/%s/%s_concatenated.a2m" \
                  % (cur_dir(), pair_id, pair_id)
out_dir = "%s/amap-multipw/%s" % (results_dir(), pair_id)
make_dir_exist(out_dir)
# Align this concatenation of the two MSAs with MUSCLE.
amap_out_path = os.path.join(out_dir, "%s_amapmulti.afa" % pair_id)
os.system("amap %s > %s" % (amap_in_path, amap_out_path))
# Extract the pairwise alignment of the seeds from the new MSA.
seedX_seq = ""
seedY_seq = ""
for record in SeqIO.parse(open(amap_out_path, "rU"), "fasta"):
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
pw_alignment_path = os.path.join(out_dir, "%s_amapmultipw.afa" % pair_id)
fp = open(pw_alignment_path, "w")
fp.write(">%s\n" % seedX_id)
fp.write("%s\n" % seedX_seq)
fp.write(">%s\n" % seedY_id)
fp.write("%s\n" % seedY_seq)
fp.close()

# If there is a structural alignment, generate CS score
reference_alignment_path = get_reference_alignment_path(pair)
if os.path.exists(reference_alignment_path):
  cs_output_file = os.path.join(out_dir, "%s_lobster_CS.out" % pair_id)
  os.system("%s -compareMSA %s -ref %s >& %s " \
	    % ( lobster_cmd(),  pw_alignment_path, reference_alignment_path, 
		cs_output_file ))
  (sp_score, cs_score) = read_sp_cs_score(cs_output_file)
  output_file = os.path.join(out_dir, "%s_amap-multipw.csv" % pair_id)
  fp = open(output_file, "w")
  fp.write("pair,SP,CS\n")
  fp.write("%s,%0.3f,%0.3f\n" % (pair_id, sp_score, cs_score))
  fp.close()

