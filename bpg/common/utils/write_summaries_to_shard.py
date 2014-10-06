#!/usr/bin/python

import os
import sys
sys.path = ['/home/ruchira/ohana_repository'] + sys.path
from bpg.common.utils.dir_of_family import get_dir_of_family_accession
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] family_accession_file"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-s", "--total_num_shards", dest="num_shards", 
                        default=1,
                        help="Number of shards into which seeds are divided")
  opt_parser.add_option("-n", "--current_shard_number", dest="shard_number",
                        default=0,
                        help="Which shard number the current job is processing")
  opt_parser.add_option("-d", "--description", dest="description",
                        default=0,
                        help="The description of download files")

  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if not os.path.exists(args[0]):
    opt_parser.error("Couldn't find %s" % args[0])
  else:
    input_file = args[0]
  try:
    num_shards = int(options.num_shards)
  except ValueError:
    opt_parser.error("--total_num_shards must be a number")
  try:
    shard_number = int(options.shard_number)
  except ValueError:
    opt_parser.error("--current_shard_number must be a number")
  f = open(input_file, "rU")
  family_accessions = [line.rstrip() for line in f.readlines()]
  f.close()
  shard_size = len(family_accessions) / (num_shards - 1)
  start = shard_number * shard_size
  end = min(start + shard_size, len(family_accessions))
  shard_file = 'downloads/%s_summaries_%03d-of-%03d' \
              % (options.description, (shard_number + 1), num_shards)
  print "Writing %s with summaries from %s to %s" % (shard_file,
                              family_accessions[start], family_accessions[end-1])
  outf = open(shard_file, 'w')
  for family_accession in family_accessions[start:end]:
      dir = get_dir_of_family_accession(family_accession)
      summary_file = os.path.join(dir, '%s_summary.dat' % family_accession)
      try:
        inf = open(summary_file)
        summary_data = inf.read()
        inf.close()
        outf.write("%s\n" % family_accession)
        outf.write(summary_data)
        outf.write("//\n")
      except IOError:
        print "Could not write %s to the shard" % summary_file
  outf.close()

if __name__ == '__main__':
  main()
