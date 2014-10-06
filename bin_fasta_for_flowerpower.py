#!/usr/bin/python

import os, sys
from optparse import OptionParser
import math
from Bio import Seq, SeqIO

def minimum_coverage(sequence_length):
  if sequence_length < 100:     
    return 0.60 + 0.0
  elif sequence_length < 200:
    return 0.60 + 0.05
  elif sequence_length < 250:
    return 0.60 + 0.10
  elif sequence_length < 300:
    return 0.60 + 0.13
  elif sequence_length < 350:
    return 0.60 + 0.15
  elif sequence_length < 400:
    return 0.60 + 0.18
  elif sequence_length < 450:
    return 0.60 + 0.20
  elif sequence_length < 500:
    return 0.60 + 0.23
  else:           # >= 500
    return 0.60 + 0.25 

def greatest_minimum_coverage():
  return 0.60 + 0.25
  
def minimum_length_that_can_cover(sequence_length):
  return int(math.ceil(minimum_coverage(sequence_length) * sequence_length))

def maximum_length_that_can_be_covered(sequence_length):
  # First check using the minimum coverage for the sequence itself
  upper_bound = math.floor(sequence_length / minimum_coverage(sequence_length))
  # However, the upper_bound is a longer sequence length, and hence may have a 
  # greater minimum coverage criterion, so sequence_length may not cover it
  # It will certainly cover a hit_length such that
  # sequence_length >= greatest_minimum_coverage() * hit_length
  lower_bound = math.floor(sequence_length / greatest_minimum_coverage())
  # Do a binary search for the greatest length which sequence_length can cover
  # This loop maintains invariant that
  # sequence_length >= minimum_coverage(lower_bound) * lower_bound
  while sequence_length < minimum_coverage(upper_bound) * upper_bound:
    middle = math.ceil((lower_bound + upper_bound)/2)
    if sequence_length >= minimum_coverage(middle) * middle:
      lower_bound = middle
    elif upper_bound <= lower_bound + 1:
      upper_bound = lower_bound
    else:
      upper_bound = middle
  return int(upper_bound)
  
def main():
  # parse command line options
  usage = "%prog [options] fasta_file_to_bin"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-i", "--bin_interval", dest="bin_interval",
              default="50",
              help="width of the range of sequence lengths to put in one bin")
  opt_parser.add_option("-m", "--max_length_to_bin", dest="max_length_to_bin",
              default="1200",
              help="length for which longer sequences all go in the last bin")
  opt_parser.add_option("-f", "--format_for_blast", dest="format_for_blast",
              action="store_true", default=False,
              help="Whether to format the binned fasta files for BLAST")
  opt_parser.add_option("-v", "--verbose", dest="verbose",
              action="store_true", default=True,
              help="Whether to print verbose output")
  opt_parser.add_option("-q", "--quiet", dest="verbose",
              action="store_false", default=True,
              help="Whether to suppress verbose output")
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  try:
    bin_interval = int(options.bin_interval)
  except ValueError:
    opt_parser.error("--bin_interval must be a number")
  try:
    max_length_to_bin = int(options.max_length_to_bin)
  except ValueError:
    opt_parser.error("--max_length_to_bin must be a number")
  global verbose
  verbose = options.verbose
  floor_num_bins, remainder = divmod(max_length_to_bin, bin_interval)
  if remainder == 0:
    num_narrow_bins = floor_num_bins
  else:
    num_narrow_bins = floor_num_bins + 1
  first_length_in_wide_bin = num_narrow_bins * bin_interval + 1
  def bin_index_of_length(sequence_length):
    if sequence_length >= first_length_in_wide_bin:
      return num_narrow_bins
    else:
      (i, rem) = divmod(sequence_length, bin_interval)
      if rem == 0:
        return i-1
      else:
        return i
  def first_length_in_bin(bin_index):
    return i * bin_interval + 1
  def last_length_in_bin(bin_index):
    return (i + 1) * bin_interval
  bins_for_hit_length = {}
  max_hit_length_binned = first_length_in_wide_bin
  file_base = os.path.splitext(os.path.split(args[0])[1])[0]
  bin_file_names = {}
  output_handles = {}
  for i in xrange(0,num_narrow_bins):
    first_seed_length = first_length_in_bin(i)
    last_seed_length = last_length_in_bin(i)
    bin_file_names[i] = "%s_%04dthrough%04d" \
      % (file_base, first_seed_length, last_seed_length)
    output_handles[i] = open(bin_file_names[i], "w")
    first_hit_length = minimum_length_that_can_cover(first_seed_length)
    last_hit_length = maximum_length_that_can_be_covered(last_seed_length)
    for j in xrange(first_hit_length,last_hit_length+1):
      if j in bins_for_hit_length.keys():
        bins_for_hit_length[j].append(i)
      else:
        bins_for_hit_length[j] = [i]
    if last_hit_length > max_hit_length_binned:
      max_hit_length_binned = last_hit_length
    if verbose:
      print "Bin #%d: %d-%d" % (i, first_seed_length, last_seed_length)
      print "Lengths that can overlap: %d-%d" \
        % (first_hit_length, last_hit_length)
  first_hit_length = minimum_length_that_can_cover(first_length_in_wide_bin)
  bin_file_names[num_narrow_bins] \
    = "%s_%4dandmore" % (file_base, first_length_in_wide_bin)
  output_handles[num_narrow_bins] = open(bin_file_names[num_narrow_bins], "w")
  if verbose:
    print "Bin #%d: %d-" % (num_narrow_bins, first_length_in_wide_bin)
    print "Lengths that can overlap: %d-" % first_hit_length
  for j in xrange(first_hit_length, max_hit_length_binned + 1):
    if j in bins_for_hit_length.keys():
      bins_for_hit_length[j].append(num_narrow_bins)
    else:
      bins_for_hit_length[j] = [num_narrow_bins]
  for record in SeqIO.parse(open(args[0], "rU"), "fasta"):
    hit_length = len(record.seq)
    if hit_length in bins_for_hit_length.keys():
      for i in bins_for_hit_length[hit_length]:
        SeqIO.write([record], output_handles[i], "fasta")
    else:
      SeqIO.write([record], output_handles[num_narrow_bins], "fasta")
  for output_handle in output_handles.values():
    output_handle.close()
  if options.format_for_blast:
    for bin_file_name in bin_file_names.values():
      cmd = 'formatdb -o T -i %s -l %s.log >& formatdb_%s.err' \
            % (bin_file_name, bin_file_name, bin_file_name)
      os.system(cmd)

if __name__ == '__main__':
  main()
