#!/usr/bin/env python

import os, commands, time
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] input_file cmd_name compute_cmd_fmt dir_cmd_fmt " \
            + "result_cmd_fmt"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-c", "--num_cores", dest="num_cores", 
                        default=1,
                        help="Number of cores to use to process inputs")
  opt_parser.add_option("-l", "--resource_list_for_qsub", dest="resource_list",
                        default="",
                        help="Resource list to pass to qsub")
  (options, args) = opt_parser.parse_args()
  if len(args) != 5:
    opt_parser.error('Incorrect number of arguments')
  if not os.path.exists(args[0]):
    opt_parser.error("Couldn't find %s" % args[0])
  else:
    input_file = args[0]
  cmd_name = args[1]
  compute_cmd_fmt = args[2]
  dir_cmd_fmt = args[3]
  result_cmd_fmt = args[4]
  try:
    num_cores = int(options.num_cores)
  except ValueError:
    opt_parser.error("--num_cores must be a number")
  f = open(input_file, "rU")
  inputs = [line.rstrip() for line in f.readlines()]
  f.close()
  cur_dir = os.getcwd()
  job_ids = {}
  input_number = 0
  for input in inputs:
    core_number = input_number % num_cores
    status, output = commands.getstatusoutput(dir_cmd_fmt % input)
    if status == 0:
      input_dir = output
      os.chdir(input_dir)
      status, output = commands.getstatusoutput(result_cmd_fmt % input)
      result_file = output.rstrip()
      if status == 0 and not os.path.exists(result_file):
        if os.path.exists("%s_%s_done" % (cmd_name, input)):
          os.system("rm %s_%s_done" % (cmd_name, input))
        input_base = os.path.split(input)[1]
        f = open("%s_%s.sh" % (cmd_name, input_base), "w")
        f.write("#!/bin/bash\n")
        f.write("source /etc/profile\n")
        f.write(compute_cmd_fmt % input_base)
        f.write("\n")
        f.write("touch %s_%s_done" % (cmd_name, input_base))
        f.close()
        cmd = "qsub -j oe -d %s -o %s_%s.out %s_%s.sh " \
                    % (input_dir, cmd_name, input_base, cmd_name, input_base)
        if options.resource_list != "":
          cmd += "-l %s " % options.resource_list
        if core_number in job_ids:
          status, output = commands.getstatusoutput("qstat %s" 
                                                      % job_ids[core_number])
          if status == 0:
            cmd += "-W depend=afterany:%s " % job_ids[core_number]
        status, output = commands.getstatusoutput(cmd)
        if status == 0:
          job_id = output.rstrip()
          f.close()
          job_ids[core_number] = job_id
      os.chdir(cur_dir)
    input_number = input_number + 1
    
if __name__ == '__main__':
  main()
