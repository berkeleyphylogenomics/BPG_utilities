#!/usr/bin/python

import sys, os
from Bio import Seq, SeqIO
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <seed_id>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]
msa_dir = "%s/gathered_homologs/%s" % (shmm_root_dir(), seed_id)
os.chdir(msa_dir)
handle = open("%s_homologs.a2m" % seed_id, "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
chain_was_cut = False
seed_key = seed_id
try:
  seed_seq = record_dict[seed_key].seq
except KeyError:
  # See if the chain identifier was cut off
  seed_key = seed_id[0:-1]
  try:
    seed_seq = record_dict[seed_key].seq
    chain_was_cut = True
  except KeyError:
    seed_key = seed_id[0:-1] + '_' + seed_id[-1:]
    seed_seq = record_dict[seed_key].seq
if chain_was_cut:
  print "The chain identifier was cut off!  Will restore..."
start_uppercase_index = 0
seeking_start_uppercase_index = True
seed_is_master = True
while (start_uppercase_index < len(seed_seq) 
        and seeking_start_uppercase_index):
  if seed_seq[start_uppercase_index].isupper():
    seeking_start_uppercase_index = False
  else:
    start_uppercase_index = start_uppercase_index + 1
  if (seed_is_master):
    if seed_seq[start_uppercase_index].islower():
      print "Warning: found lowercase character " + \
      "before any uppercase character in aligned seed sequence"
      print "The seed is not the master in this alignment"
      seed_is_master = False
if (seeking_start_uppercase_index):
  print "Unable to find any uppercase characters in aligned seed sequence"
  print "Exiting..."
  sys.exit(0)
end_uppercase_index = len(seed_seq) - 1
seeking_end_uppercase_index = True
while (end_uppercase_index >= 0 and seeking_end_uppercase_index):
  if seed_seq[end_uppercase_index].isupper():
    seeking_end_uppercase_index = False
  else:
    end_uppercase_index = end_uppercase_index - 1
  if (seed_is_master):
    if seed_seq[end_uppercase_index].islower():
      print "Warning: found lowercase character " + \
      "after all uppercase characters in aligned seed sequence"
      print "The seed is not the master in this alignment"
      seed_is_master = False
# want to use the index as the end of a range
end_uppercase_index = end_uppercase_index + 1
cropped_homolog_dir = os.path.join(shmm_root_dir(), "cropped_gathered_homologs")
output_dir = os.path.join(cropped_homolog_dir, seed_id)
make_dir_exist(output_dir)
output_file_name = os.path.join(output_dir, 
                                "%s_homologs_cropped_to_seed.a2m" % seed_id)
print "start_uppercase_index: %d" % start_uppercase_index
print "end_uppercase_index: %d" % end_uppercase_index
fp = open(output_file_name, "w")
for key in record_dict.keys():
  if key == seed_key:
    fp.write(">%s\n" % seed_id)
  else:
    fp.write(">%s\n" % key)
  seq = record_dict[key].seq
  fp.write("%s\n" % seq[start_uppercase_index:end_uppercase_index].tostring())
fp.close()
print "Wrote %s" % output_file_name
