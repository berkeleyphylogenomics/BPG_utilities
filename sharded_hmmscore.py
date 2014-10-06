#!/usr/bin/env python

import os
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] sharded_fasta_file_basepath hmm_path"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-n", "--num_shards", dest="num_shards",
        default="20", help="number of shards into which to split fasta file")
  (options, args) = opt_parser.parse_args()
  if len(args) != 2:
    opt_parser.error("Incorrect number of arguments")
  try:
    num_shards = int(options.num_shards)
  except ValueError:
    opt_parser.error("--num_shards must be a positive number")
  if num_shards <= 0:
    opt_parser.error("--num_shards must be a positive number")
  base_path = args[0]
  base_name = os.path.split(base_path)[1]
  hmm_path = args[1]
  hmm_base = os.path.splitext(os.path.split(hmm_path)[1])[0]
  if not os.path.exists(hmm_path):
    opt_parser.error("HMM file %s does not exist" % hmm_path)
  shard_filepaths = {}
  for shard_num in xrange(0, num_shards):
    shard_filepaths[shard_num] \
      = "%s_%d-of-%d" % (base_path, shard_num, num_shards)
    if not os.path.exists(shard_filepaths[shard_num]):
      opt_parser.error("Fasta shard file %s does not exist"
                        % shard_filepaths[shard_num])
  for shard_num in xrange(0, num_shards):
    runname = "%s_%d-of-%d_vs_%s" % (base_name, shard_num, num_shards,
                                      hmm_base)
    cmd = "hmmscore %s -dbsize 100000 -db %s -sw 2 -adpstyle 5 -i %s " \
          % (runname, shard_filepaths[shard_num], hmm_path)
    cmd = cmd + "> hmmscore.%s.out\n" % runname
    job_file = "do_hmmscore_%s.sh" % runname
    f = open(job_file, "w")
    f.write("#!/bin/bash\n")
    f.write(cmd)
    f.write("touch %s_is_done\n" % runname)
    f.close()
    os.system("qsub -j oe -d %s -o do_hmmscore_%s.out %s" 
              % (os.getcwd(), runname, job_file))

if __name__ == '__main__':
  main()
