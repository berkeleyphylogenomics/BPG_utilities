#!/usr/bin/env python

import commands
import glob
from optparse import OptionParser
import os
import string
import time

from bpg.common.submit_job_on_queue import submit_job_on_queue

def main():
  # parse command line options
  usage = "%prog [options] sharded_fasta_file_basepath hmm_path"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-n", "--num_shards", dest="num_shards",
        default="20", help="number of shards into which to split fasta file")
  opt_parser.add_option("-e", "--max_evalue", dest="max_evalue",
        default="0.001", help="maximum e-value to include in final score file")
  (options, args) = opt_parser.parse_args()
  if len(args) != 2:
    opt_parser.error("Incorrect number of arguments")
  try:
    num_shards = int(options.num_shards)
  except ValueError:
    opt_parser.error("--num_shards must be a positive number")
  if num_shards <= 0:
    opt_parser.error("--num_shards must be a positive number")
  try:
    max_evalue = float(options.max_evalue)
  except ValueError:
    opt_parser.error("--max_evalue must be a nonnegative real")
  if max_evalue < 0:
    opt_parser.error("--max_evalue must be a nonnegative real")
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
  hmmscore_cmd_of_shard = {}
  job_id_of_hmmscore = {}
  shards_to_score = set()
  def runname_of_shard(shard_num):
    return "%s_%d-of-%d_vs_%s" % (base_name, shard_num, num_shards, hmm_base)
  for shard_num in xrange(0, num_shards):
    runname = runname_of_shard(shard_num)
    shards_to_score.add(runname)
    cmd = "hmmscore %s -dbsize 100000 -db %s -sw 2 -adpstyle 5 -i %s " \
          % (runname, shard_filepaths[shard_num], hmm_path)
    cmd = cmd + "> hmmscore.%s.out\n" % runname
    hmmscore_cmd_of_shard[runname] = cmd
    job_id_of_hmmscore[runname] = submit_job_on_queue('hmmscore', cmd, runname)

  time.sleep(1)

  while len(shards_to_score) > 0:
    done_shards = glob.glob("hmmscore_*_done")
    for done_file in done_shards:
      shards_to_score.remove(done_file[9:-5])
      os.system("rm %s" % done_file)
    for runname in shards_to_score:
      status, output = commands.getstatusoutput("qstat %s" %
                                                job_id_of_hmmscore[runname])
      if status != 0 and not os.path.exists("hmmscore_%s_done" % runname):
        job_id_of_hmmscore[runname] = submit_job_on_queue('hmmscore',
                                      hmmscore_cmd_of_shard[runname], runname)

  of = open("%s_vs_%s.dist" % (base_name, hmm_base), "w")
  first_distfile = True
  dist_lines = []
  for runname in hmmscore_cmd_of_shard.keys():
    distfile = "%s.dist" % runname
    f = open(distfile)
    lines = f.readlines()
    f.close()
    for line in lines:
      if line == "":
        break
      if line[0] == '%':
        if first_distfile:
          if line[0:17] == '% Database Files:':
            of.write('%% Database Files:  %s\n' % base_path)
          else:
            of.write(line)
      else:
        e_value = string.atof(line.split()[4])
        dist_lines.append((e_value, line))
    first_distfile = False
  dist_lines.sort()
  for e_value, line in dist_lines:
    if e_value > max_evalue:
      break
    of.write(line)
  of.close()

if __name__ == '__main__':
  main()
