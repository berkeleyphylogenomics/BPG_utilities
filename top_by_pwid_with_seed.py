#!/usr/bin/python

import sys, os, string
from Bio import Seq, SeqIO
import BPG_common.BPGPWID
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 3:
  print "Usage: %s <num_to_keep> <seed_id>" % sys.argv[0]
  sys.exit(0)

num_to_keep = int(sys.argv[1])
seed_id = sys.argv[2]
shmmshmmroot = shmm_root_dir()
msa_dir = "%s/nr100_cropped_gathered_homologs/%s" % (shmmshmmroot, seed_id)
os.chdir(msa_dir)
new_dir = "%s/top%d_nr100_cropped_gathered_homologs/%s" \
          % (shmmshmmroot, num_to_keep, seed_id)
make_dir_exist(new_dir)
# We need to look in the cropped msa prior to unique90 for the seed sequence,
# as unique90 may have dropped it
cropped_msa_dir = "%s/cropped_gathered_homologs/%s" % (shmmshmmroot, seed_id)
seed_seq = ""
handle = open(os.path.join( cropped_msa_dir, 
                            "%s_homologs_cropped_to_seed.a2m" % seed_id), "rU")
for record in SeqIO.parse(handle, "fasta"):
  if record.id == seed_id or record.id == seed_id[0:-1]:
    seed_seq = record.seq
handle.close()
if seed_seq == "":
  print "Couldn't find the seed sequence, exiting"
  sys.exit(0)
seed_str = seed_seq.tostring()
handle = open("%s_homologs_cropped_nr100.a2m" % seed_id, "rU")
record_list = list(SeqIO.parse(handle, "fasta"))
handle.close()
pwid_dict = {}
for record in record_list:
  pwid_dict[record.id] = \
    BPG_common.BPGPWID.pairwise_identity_KS_1(seed_str, record.seq.tostring())
def compare_seqs(record0, record1):
  return cmp(pwid_dict[record0.id], pwid_dict[record1.id])
record_list.sort(compare_seqs, reverse=True)
fp = open(os.path.join(new_dir, "%s.a2m" % seed_id), "w")
fp.write(">%s\n" % seed_id)
fp.write("%s\n" % seed_str)
i = 0
for record in record_list:
  if i < num_to_keep:
    if record.id != seed_id and record.id != seed_id[0:-1]:
      fp.write(">%s\n" % record.id)
      fp.write("%s\n" % record.seq.tostring())
    i = i + 1
fp.close()
print "Wrote %s/%s.a2m" % (new_dir, seed_id)
