#!/usr/bin/python

import os, sys, string
from Bio import Seq, SeqIO
from matchmaker.shmm_shmm_lib import *

seedX_id = sys.argv[1]
seedY_id = sys.argv[2]
pair_id = "%s_%s" % (seedX_id, seedY_id)

# FlowerPower has gathered homologs for each seed and created an alignment.
# ClustalW complains if there is more than one sequence with the same name.
# Concatenate the two alignments, renaming to avoid duplicate headers
seedX_alignment_path = "%s/%s/%s.a2m" % (cur_dir(), seedX_id, seedX_id)
seedY_alignment_path = "%s/%s/%s.a2m" % (cur_dir(), seedY_id, seedY_id)
recordX_dict = SeqIO.to_dict(SeqIO.parse(
                                    open(seedX_alignment_path, "r"), "fasta"))
recordY_dict = SeqIO.to_dict(SeqIO.parse(
                                    open(seedY_alignment_path, "r"), "fasta"))
for headerX in recordX_dict.keys():
  if headerX in recordY_dict:
    headerY_components = recordY_dict[headerX].id.split('|')
    headerY_components.insert(1, "AlignedTo%s" % seedY_id)
    recordY_dict[headerX].id = '|'.join(headerY_components)
for headerY in recordY_dict.keys():
  if headerY in recordX_dict:
    headerX_components = recordX_dict[headerY].id.split('|')
    headerX_components.insert(1, "AlignedTo%s" % seedX_id)
    recordX_dict[headerY].id = '|'.join(headerX_components)
out_dir = "%s/clustalw-multipw/%s" % (shmm_root_dir(), pair_id)
make_dir_exist(out_dir)
clustalw_in_path = os.path.join(out_dir, "%s_concat_for_clustalw.a2m" % pair_id)
fp = open(clustalw_in_path, "w")
for key in recordX_dict.keys():
  fp.write(">%s\n" % recordX_dict[key].id)
  fp.write("%s\n" % recordX_dict[key].seq.tostring())
for key in recordY_dict.keys():
  fp.write(">%s\n" % recordY_dict[key].id)
  fp.write("%s\n" % recordY_dict[key].seq.tostring())
fp.close()
# Align this concatenation of the two MSAs with ClustalW.
clustalw_out_path = os.path.join(out_dir, "%s_clustalwmulti.gde" % pair_id)
cmd = "clustalw -infile=%s -outfile=%s -output=GDE" \
          % (clustalw_in_path, clustalw_out_path)
#print "Running", cmd
os.system(cmd)

# Make GDE file into aligned FASTA file
msa_alignment_filename = "%s_clustalwmulti.msa" % pair_id
msa_file_path = os.path.join(out_dir, msa_alignment_filename)
ifp = open(clustalw_out_path, "r")
lines = ifp.readlines()
ifp.close()
ofp = open(msa_file_path, "w")
for line in lines:
  if line[0] == '%':
    ofp.write(string.replace(line,"%",">"))
  else:
    ofp.write(string.upper(line))
ofp.close()
# Extract the pairwise alignment of the seeds from the new MSA.
seedX_seq = ""
seedY_seq = ""
for record in SeqIO.parse(open(msa_file_path, "rU"), "fasta"):
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
pw_alignment_filename = "%s_clustalwmultipw.afa" % pair_id
pw_alignment_path = os.path.join(out_dir, pw_alignment_filename)
"""
fp = open(pw_alignment_path, "w")
fp.write(">%s\n" % seedX_id)
fp.write("%s\n" % seedX_seq)
fp.write(">%s\n" % seedY_id)
fp.write("%s\n" % seedY_seq)
fp.close()
"""
write_alignment_file_compressed(seedX_seq, seedY_seq, pw_alignment_path)

# If there is a structural alignment, generate the CS score.
reference_alignment_path = os.path.join(pair_dir(),
                                        pair_id, "%s.afa" % pair_id)
if os.path.exists(reference_alignment_path):
  cs_output_file = os.path.join(out_dir, "%s_lobster_CS.out" % pair_id)
  os.system("%s -compareMSA %s -ref %s >& %s " \
	    % ( lobster_cmd(), pw_alignment_path, reference_alignment_path, 
		cs_output_file ))
  (sp_score, cs_score) = read_sp_cs_score(cs_output_file)

  output_file = os.path.join(out_dir, "%s_clustalw-multipw.csv" % pair_id)
  fp = open(output_file, "w")
  fp.write("pair,SP,CS\n")
  fp.write("%s,%0.3f,%0.3f\n" % (pair_id, sp_score, cs_score))
  fp.close()

