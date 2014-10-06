#!/usr/bin/python

import os, glob
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s  <seed_id>" % sys.argv[0]
  sys.exit(0)

id = sys.argv[1]
shmmshmmroot = shmm_root_dir()
os.chdir(shmmshmmroot)
make_dir_exist("gathered_homologs/%s" % id)
os.chdir("gathered_homologs/%s" % id)
if not os.path.exists("%s_homologs.a2m" % id):
  cmd = "ln -s ../../flowerpower_seedmaster/%s/final.a2m %s_homologs.a2m" \
          % (id, id)
  os.system(cmd)
