#!/usr/bin/python
import os, sys, glob, cPickle, time
import subprocess, shlex
from optparse import OptionParser
from Bio import Seq, SeqIO

def main():
  # parse command line options
  usage = "%prog [options] input_file runname compute_cmd_fmt dir_cmd_fmt " \
            + "result_cmd_fmt"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-s", "--total_num_shards", dest="num_shards", 
                        default=1,
                        help="Number of shards into which seeds are divided")
  opt_parser.add_option("-n", "--current_shard_number", dest="shard_number",
                        default=0,
                        help="Which shard number the current job is processing")
  (options, args) = opt_parser.parse_args()
  if len(args) != 5:
    opt_parser.error('Incorrect number of arguments')
  if not os.path.exists(args[0]):
    opt_parser.error("Couldn't find %s" % args[0])
  else:
    input_file = args[0]
  runname = args[1]
  compute_cmd_fmt = args[2]
  dir_cmd_fmt = args[3]
  result_cmd_fmt = args[4]
  try:
    num_shards = int(options.num_shards)
  except ValueError:
    opt_parser.error("--total_num_shards must be a number")
  try:
    shard_number = int(options.shard_number)
  except ValueError:
    opt_parser.error("--current_shard_number must be a number")
  cur_dir = os.getcwd()
  f = open(input_file, "rU")
  inputs = [line.rstrip() for line in f.readlines()]
  f.close()
  input_number = 0
  for input in inputs:
    if input_number % num_shards == shard_number:
      pipe = subprocess.Popen(dir_cmd_fmt % input, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
      (output, errout) = pipe.communicate(input=None)
      status = pipe.returncode
      if status == 0:
        input_dir = output.rstrip()
        os.chdir(input_dir)
        pipe = subprocess.Popen(result_cmd_fmt % input, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, errout) = pipe.communicate(input=None)
        status = pipe.returncode
        result_file = output.rstrip()
        if status == 0 and not os.path.exists(result_file):
          if os.path.exists("%s_%s_done" % (runname, input)):
            subprocess.call(shlex.split("rm %s_%s_done" % (runname, input)))
          input_base = os.path.split(input)[1]
          cmd = compute_cmd_fmt % input_base
          sys.stdout.write("Running on input #%d: %s$ %s\n" 
                            % (input_number, input_dir, cmd))
          sys.stdout.flush()
          status = subprocess.call(shlex.split(cmd))
          if status != 0:
            sys.stdout.write("Nonzero exit status %d\n" % status)
            sys.stdout.flush()
          else:
            sys.stdout.write("Normal exit\n")
            sys.stdout.flush()
          subprocess.call(shlex.split("touch %s_%s_done" 
                                      % (runname, input_base)))
          if not os.path.exists(result_file):
            sys.stdout.write("%s for %s FAILED!\n" % (runname, input_base))
            sys.stdout.flush()
        elif status != 0:
          print "%s exited with nonzero status %d" % (result_cmd_fmt % input,
                                                      status)
          print result_file
        os.chdir(cur_dir)
    input_number = input_number + 1
  subprocess.call(shlex.split('touch %s_job_%d_of_%d_done' 
                              % (runname, shard_number, num_shards)))

if __name__ == '__main__':
  main()
