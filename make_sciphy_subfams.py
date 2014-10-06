#!/usr/bin/python

# Make SCI-PHY subfamilies and SHMMSs from the MSA for this seed.

import os, sys
from matchmaker.shmm_shmm_lib import *

out_filename = "shmm.out"

if len(sys.argv) < 2:
  print "Usage: %s <seed_id>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]
root = "%s/%s" % (cur_dir(), seed_id)
alignment_path = os.path.join(root, "%s.a2m" % seed_id)
sciphy_dir = os.path.join(root, "sciphy")
make_dir_exist(sciphy_dir)
os.chdir(sciphy_dir)

if sciphy_hmms_completed(seed_id, sciphy_dir, out_filename):
  print "SCI-PHY subfams/hmms previously created for %s" % seed_id

else:
  # Link MSA and ascii GHMM to this directory
  os.system("ln -s ../%s.a2m %s.a2m" % (seed_id, seed_id))
  os.system("ln -s ../ascii_ghmm/%s.mod %s.mod" % (seed_id, seed_id))

  # Run SCI-PHY to create subfamilies and SHMMs
  cmd = "%s %s -i %s.a2m -hmm %s.mod >& %s" % (sciphy_cmd(), seed_id,
       seed_id, seed_id, out_filename)
  print cmd
  os.system( cmd )
 
  # Rename the resulting SHMMs
  cmd = 'mmv -d "%s.N*.mod" "%s.sciphy.sf#1.mod"' % (seed_id, seed_id)
  print cmd
  os.system( cmd )
