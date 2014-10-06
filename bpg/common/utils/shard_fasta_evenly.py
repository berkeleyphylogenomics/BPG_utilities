#!/usr/bin/env python

import os, sys
from stat import *
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] fasta_file_to_shard"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-n", "--num_shards", dest="num_shards",
        default="20", help="number of shards into which to split fasta file")
  opt_parser.add_option("-b", "--make_blastable", dest="make_blastable",
        action="store_true", default=False,
        help="whether to make each shard into a blastable database")
  opt_parser.add_option("--nomake_blastable", dest="make_blastable",
        action="store_false", default=False,
        help="whether not to make each shard into a blastable database")
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  try:
    num_shards = int(options.num_shards)
  except ValueError:
    opt_parser.error("--num_shards must be a positive number")
  if num_shards <= 0:
    opt_parser.error("--num_shards must be a positive number")
  fasta_path = args[0]
  file_base = os.path.splitext(os.path.split(fasta_path)[1])[0]
  file_size = os.stat(fasta_path)[ST_SIZE]
  f = open(fasta_path, "r")
  target_shard_size = file_size / num_shards
  eof = False
  for shard_num in xrange(0, num_shards):
    shard_filename = "%s_%d-of-%d" % (file_base, shard_num, num_shards)
    out_file = open(shard_filename, "w")
    out_file.write(f.read(target_shard_size))
    out_file.write(f.readline())
    seeking_next_header = True
    while seeking_next_header:
      start_of_line_pos = f.tell()
      line = f.readline()
      if len(line) > 0:
        if line[0] == '>':
          seeking_next_header = False
          f.seek(start_of_line_pos)
        else:
          out_file.write(line)
      else:
        seeking_next_header = False
        eof = True
    out_file.close()
    if options.make_blastable:
      cmd = 'formatdb -o F -i %s -l %s.log >& %s.err' \
          % (shard_filename, shard_filename, shard_filename)
      os.system(cmd)
    if eof:
      if shard_num < num_shards - 1:
        print "Sharding terminated early: total %d shards" % (shard_num + 1)
      break
  f.close()

  

if __name__ == '__main__':
  main()
