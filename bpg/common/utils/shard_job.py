#!/usr/bin/python

import os, sys, string

def usage():
  print "Usage: %s <runname> <cmd> <num_shards> [arguments to qsub]" % sys.argv[0]
if len(sys.argv) < 4:
  usage()
  sys.exit(0)

try:
  num_shards = string.atoi(sys.argv[3])
except ValueError:
  usage()
  sys.exit(0)

runname = sys.argv[1]
cmd = sys.argv[2]

num_digits = "%d" % len("%d" % num_shards)

for i in xrange(0,num_shards):
  filebase_fmt = "%s_job_%0" + num_digits + "d_of_%0" + num_digits + "d"
  filebase = filebase_fmt % (runname, i, num_shards)
  f = open("%s.sh" % filebase, "w")
  f.write("#!/bin/bash\n")
  f.write("source /clusterfs/ohana/software/ohana_environment.sh\n")
  f.write("development\n")
  f.write("%s -n %d -s %d\n" % (cmd, i, num_shards))
  f.close()
  qsubcmd = "qsub -j oe -d %s -o %s.out %s -lwalltime=288:00:00 %s.sh" \
            % (os.getcwd(), filebase, ' '.join(sys.argv[4:]), filebase)
  print qsubcmd
  os.system(qsubcmd)
