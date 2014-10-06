#!/usr/bin/python

import os, sys
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 3:
  print "Usage: %s <seedX_id> <seedY_id>" % sys.argv[0]
  sys.exit(0)

seedX_id = sys.argv[1]
seedY_id = sys.argv[2]
pair_id = "%s_%s" % (seedX_id, seedY_id)

# Concatenate the alignments of the homologs for each seed
root = cur_dir()
seedX_alignment_path = "%s/%s/%s.a2m" % (root, seedX_id, seedX_id)
seedY_alignment_path = "%s/%s/%s.a2m" % (root, seedY_id, seedY_id)
out_dir = "%s/concat_for_multipw/%s" % (root, pair_id)
make_dir_exist(out_dir)
out_path = os.path.join(out_dir, "%s_concatenated.a2m" % pair_id)
os.system("cat %s %s > %s" % (seedX_alignment_path, seedY_alignment_path,
                              out_path))
