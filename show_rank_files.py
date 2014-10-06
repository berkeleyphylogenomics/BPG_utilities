#!/bin/env python

"""
Show files created by uniquification, clustering, and ranking
for a target/template pair.
"""

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <pair>" % sys.argv[0]
  sys.exit(1)

(seedX, seedY) = sys.argv[1].split("_")
dir = align_dir(seedX, seedY)
filename = "%s_%s*" % (seedX, seedY)
filepath = os.path.join(dir, filename)
cmd = "ls -rltd %s" % (filepath)
os.system(cmd)
cmd = "wc %s" % (filepath)
os.system(cmd)
