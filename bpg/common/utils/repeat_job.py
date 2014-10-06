#!/usr/bin/python

import os, sys, string

def usage():
  print "Usage: %s <runname> <cmd> <num_shards> [arguments for qsub]" % sys.argv[0]
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

f = open("%s.sh" % runname, "w")
f.write("#!/bin/bash\n")
f.write("source /clusterfs/ohana/software/ohana_environment.sh\n")
f.write("development\n")
f.write(cmd)
f.write("\n")
f.close()
os.system("chmod a+x %s.sh" % runname)
for i in xrange(0,num_shards):
  qsubcmd = "qsub -S /bin/bash -j oe -d %s %s.sh %s" \
            % (os.getcwd(), runname, ' '.join(sys.argv[4:]))
  print qsubcmd
  os.system(qsubcmd)
