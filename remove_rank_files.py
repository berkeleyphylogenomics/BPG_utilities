#!/bin/env python

"""
Remove all files created by Sriram's ranking method
for a target/template pair.
"""

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 3:
  print "Usage: %s <template_id> <target_id>>" % sys.argv[0]
  sys.exit(1)

seedX = sys.argv[1]
seedY = sys.argv[2]
dir = align_dir(seedX, seedY)
suffix_list = [
	".rank",
	".ss2",
	".yl",
	"_ss2.csv",
	".features",
	".features.ext",
]

for suffix in suffix_list:
  filename = "%s_%s%s" % (seedX, seedY, suffix)
  filepath = os.path.join(dir, filename)
  cmd = "rm %s" % (filepath)
  os.system(cmd)
