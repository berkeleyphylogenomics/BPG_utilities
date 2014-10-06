#!/usr/bin/python

import os, sys
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <seed_id>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]
root = "%s/%s" % (cur_dir(), seed_id)
os.chdir(root)
os.system("w0.5 %s.a2m %s.mod" % (seed_id, seed_id))
