#!/usr/bin/python

import sys, os
from Bio import Seq, SeqIO
import BPG_common.BPGPWID
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <query_id>" % sys.argv[0]
  sys.exit(0)

query_id = sys.argv[1]
msa_dir = "%s/flowerpower_seedmaster/%s" % (shmm_root_dir(), query_id)
os.chdir(msa_dir)
new_dir = "%s/uniqueseq_top1000/%s" % (shmm_root_dir(), query_id)
make_dir_exist(new_dir)
query_seq = ""
for record in SeqIO.parse(open("final.a2m", "rU"), "fasta"):
  if record.id == query_id or record.id == query_id[0:-1]:
    query_seq = record.seq
if query_seq == "":
  print "Couldn't find the query sequence, exiting"
  sys.exit(0)
query_str = query_seq.tostring()
handle = open("%s_unique90.a2m" % query_id, "rU")
record_list = list(SeqIO.parse(handle, "fasta"))
pwid_dict = {}
for record in record_list:
  pwid_dict[record.id] = \
    BPG_common.BPGPWID.pairwise_identity_KS_1(query_str, record.seq.tostring())
def compare_seqs(record0, record1):
  return cmp(pwid_dict[record0.id], pwid_dict[record1.id])
record_list.sort(compare_seqs)
fp = open(os.path.join(new_dir, "%s_uniqueseq_top1000.a2m" % query_id), "w")
fp.write(">%s\n" % query_id)
fp.write("%s\n" % query_str)
i = 0
for record in record_list:
  if i < 1000:
    if record.id != query_id and record.id != query_id[0:-1]:
      fp.write(">%s\n" % record.id)
      fp.write("%s\n" % record.seq.tostring())
      i = i + 1
fp.close()
