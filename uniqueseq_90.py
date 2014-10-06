#!/usr/bin/python
import sys, os
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <seed_id>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]

indir = "%s/cropped_gathered_homologs/%s" % (shmm_root_dir(), seed_id)
          #% seed_id
outdir = "%s/unique90_cropped_gathered_homologs/%s"  % (shmm_root_dir(),
     seed_id)
make_dir_exist(outdir)
os.chdir(outdir)
os.system("uniqueseq %s_cropped_unique90 -percentid 0.90 " % seed_id + \
    "-alignfile %s/%s_homologs_cropped_to_seed.a2m"  % (indir, seed_id))
