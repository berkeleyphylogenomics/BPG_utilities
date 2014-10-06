#!/usr/bin/env python

# # This version of top_by_pwid is called by flowerpower.pl
# # to trim excess homologs when appropriate.


import sys, os, string
from optparse import OptionParser
from Bio import Seq, SeqIO
import BPG_common.BPGPWID


def main():

  parser = OptionParser()
  (options, args) = parser.parse_args()

  if len(args) != 4:
    print "Usage: %s <input file (FASTA format)> <output file> <number to keep> <seed ID>" % sys.argv[0]
    sys.exit(1)

  input_filepath  = sys.argv[1]
  output_filepath = sys.argv[2]
  try:
    num_to_keep = int(sys.argv[3])
  except ValueError:
    print "The value [%s] appears to not be an integer.  Exiting..." % (sys.argv[3])
    sys.exit(1)
  seed_id = sys.argv[4]
  
  print input_filepath
  print output_filepath
  print num_to_keep
  print seed_id

  seed_seq = ""

  # # Establish I/O
  if not os.path.exists(input_filepath) :
    print "Couldn't find the input file [%s].  Exiting..." % (input_filepath)
    sys.exit(1)

  if not os.access(input_filepath, os.R_OK):
    print "Couldn't read the input file [%s].  Exiting..." % (input_filepath)
    sys.exit(1)
      
  input_handle = open(input_filepath, "rU")

  output_dir_and_file = os.path.split(output_filepath)
  output_dir  = output_dir_and_file[0]
  output_file = output_dir_and_file[1]

  if "" == output_dir :
    output_dir = "./"

  if not os.access(output_dir, os.W_OK):
    print "No access to write output in directory [%s].  Exiting..." % (output_dir)
    sys.exit(1)

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
    sys.exit(2)

  # Remove the dots from the seed sequence, since uniqueseq has removed the dots
  # from everything else
  # seed_str = string.replace(seed_seq.tostring(), '.', '')
  seed_str = seed_seq

  # # Key: Record ID
  # # Value: % identity with the seed
  pwid_dict = {}

  # # Enter each input sequence in the dictionary.
  for record in record_list:
    pwid_dict[record.id] = \
      BPG_common.BPGPWID.pairwise_identity_KS_1(seed_str, record.seq.tostring())

  # # For a pair of sequences, compare their % identities with the seed (for sorting)
  def compare_seqs(record0, record1):
    return cmp(pwid_dict[record0.id], pwid_dict[record1.id])

  # # Sort the list of sequences by % identity with the seed.
  record_list.sort(compare_seqs, reverse=True)

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

  # print "Wrote %s" % (output_filepath)
  sys.exit(0)

if __name__=="__main__":
  main()
