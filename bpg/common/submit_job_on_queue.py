#!/usr/bin/env python

import os
import commands

def submit_job_on_queue(progstr, cmd, input):
  f = open("%s_%s.sh" % (progstr, input), "w")
  f.write("#!/bin/bash\n")
  f.write("touch %s_%s_start\n" % (progstr, input))
  f.write("%s\n" % cmd)
  f.write("touch %s_%s_done\n" % (progstr, input))
  f.close()
  qsubcmd = 'ssh -v ohana.berkeley.edu "qsub -j oe ' + \
            '-d %s -o %s_%s.out %s_%s.sh ' \
            % (os.getcwd(), progstr, input, progstr, input) + \
            '>& %s/qsub_%s_%s.out"' % (os.getcwd(), progstr, input)
  for i in xrange(0,3):
    status, output = commands.getstatusoutput(qsubcmd)
    if status == 0:
      f = open("%s/qsub_%s_%s.out" % (os.getcwd(), progstr, input))
      job_id = f.readline().rstrip()
      f.close()
      break
    else:
      print "ssh failed"
      print output
  if status != 0:
    raise "%s '%s' failed" % (progstr, qsubcmd)
  return job_id
