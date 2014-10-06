#!/usr/bin/python

import os
import string
import sys

from Bio import Seq, SeqIO

from bpg.common import BPGPWID

#from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 5:
  print "Usage: %s <input file (FASTA format)> <outut file> <number to keep> <seed id>" % sys.argv[0]
  sys.exit(0)

input_filepath = sys.argv[1]
output_filepath = sys.argv[2]

# if not isinstance(sys.argv[3], int):
#   print "The value [%s] appears to not be an integer.  Exiting..." % (sys.argv[3])
#   sys.exit(0)

try:
  num_to_keep = int(sys.argv[3])
except:
  print "The value [%s] appears to not be an integer.  Exiting..." % (sys.argv[3])
  sys.exit(0)
  

seed_id = sys.argv[4]

#shmmshmmroot = shmm_root_dir()
#msa_dir = "%s/nr100_cropped_gathered_homologs/%s" % (shmmshmmroot, seed_id)
#os.chdir(msa_dir)
#new_dir = "%s/top%d_nr100_cropped_gathered_homologs/%s" \
#          % (shmmshmmroot, num_to_keep, seed_id)
#make_dir_exist(new_dir)
# We need to look in the cropped msa prior to unique90 for the seed sequence,
# as unique90 may have dropped it
#cropped_msa_dir = "%s/cropped_gathered_homologs/%s" % (shmmshmmroot, seed_id)

seed_seq = ""

#handle = open(os.path.join( cropped_msa_dir, 
#                            "%s_homologs_cropped_to_seed.a2m" % seed_id), "rU")

# # Establish I/O
if not os.path.exists(input_filepath) :
  print "Couldn't find the input file [%s].  Exiting..." % (input_filepath)
  sys.exit(0)

if not os.access(input_filepath, os.R_OK):
  print "Couldn't read the input file [%s].  Exiting..." % (input_filepath)
  sys.exit(0)

input_handle = open(input_filepath, "rU")

output_dir_and_file = os.path.split(output_filepath)
output_dir  = output_dir_and_file[0]
output_file = output_dir_and_file[1]

if "" == output_dir :
  output_dir = "./"

if not os.access(output_dir, os.W_OK):
  print "No access to write output in directory [%s].  Exiting..." % (output_dir)
  sys.exit(0)

# # Read and parse the input
record_list = list(SeqIO.parse(input_handle, "fasta"))

input_handle.close()

# # Scan for the seed sequence
for record in record_list:
  if record.id == seed_id or record.id == seed_id[0:-1]:
    seed_seq = record.seq
    break

if seed_seq == "":
  print "Couldn't find the seed sequence.  Exiting..."
  sys.exit(0)

# Remove the dots from the seed sequence, since uniqueseq has removed the dots
# from everything else
seed_str = string.replace(seed_seq.tostring(), '.', '')

# # Key: Record ID
# # Value: % identity with the seed
pwid_dict = {}

# # Enter each input sequence in the dictionary.
for record in record_list:
  pwid_dict[record.id] = \
    BPGPWID.pairwise_identity_KS_1(seed_str, record.seq.tostring())

# # For a pair od sequences, compare their % identities with the seed (for sorting)
def compare_seqs(record0, record1):
  return cmp(pwid_dict[record0.id], pwid_dict[record1.id])

# # Sort the list of sequences by % identity with the seed.
record_list.sort(compare_seqs)

output_handle = open(output_filepath, "wU")
output_handle.write(">%s\n" % seed_id)
output_handle.write("%s\n" % seed_str)

# # Write the output.
num_written = 1
for record in record_list:
  if num_written < num_to_keep:
    if record.id != seed_id and record.id != seed_id[0:-1]:
      output_handle.write(">%s\n" % record.id)
      output_handle.write("%s\n" % record.seq.tostring())
      num_written = num_written + 1

output_handle.close()

print "Wrote %s" % (output_filepath)
