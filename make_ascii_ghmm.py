#!/usr/bin/python

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <seed_id>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]
main_dir = os.path.join(cur_dir(), seed_id)
os.chdir(main_dir)
ascii_ghmm_dir = os.path.join(main_dir, "ascii_ghmm")
make_dir_exist(ascii_ghmm_dir)
os.system("cp %s.mod %s" % (seed_id, ascii_ghmm_dir))
os.chdir(ascii_ghmm_dir)
os.system("hmmconvert %s -model_file %s.mod > hmmconvert.out" \
          % (seed_id, seed_id))
os.chdir(main_dir)
kerfdirs = glob.glob("kerf[1-9][0-9]")
for kerfdir in kerfdirs:
  kerfsharedir = os.path.join(kerfdir, "info-share")
  make_dir_exist(kerfsharedir)
  os.chdir(kerfsharedir)
  if not os.path.exists("%s.mod" % seed_id):
    os.symlink("../../ascii_ghmm/%s.mod" % seed_id, "%s.mod" % seed_id)
  os.chdir(main_dir)
