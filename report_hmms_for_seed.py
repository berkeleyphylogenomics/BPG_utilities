#!/usr/bin/python

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

"""
report_hmms_for_seed.py  June 3, 2008  Terry Farrah  Sjolander Lab

Part of the matchmaker (SHMM-SHMM) alignment pipeline.
Counts the SHMMs found in the directory tree for a seed sequence.
Only counts SHMMs in directories hard-coded below.
Used to evaluate, at the end of the single-seq phase of the pipeline,
whether the pipeline completed successfully.

Can call as function. Returns True if any HMMs are found in
the proper place.
"""

def report_hmms_for_seed(seed_id, verbose=False):

  initial_wd = os.getcwd()
  shmm_spec = "%s.*.sf*.mod" % (seed_id)
  seed_dir = os.path.join(cur_dir(), seed_id)
  if os.path.exists(seed_dir):
    os.chdir(seed_dir)
  else:
    if verbose: print >> sys.stderr, \
         "%s: directory %s does not exist." % (sys.argv[0], seed_dir)
    return(False)

  hmms_found = False

  # Report kerf cuts
  kerf_dirs = glob.glob("kerf*")
  for kerf_dir in kerf_dirs:
    percent_id = int(kerf_dir[4:])
    os.chdir(kerf_dir)
    num_subfams  = len(glob.glob("*.fa"))
    if num_subfams > 0:  hmms_found = True
    if verbose: print "%s: %d subfamilies" % (kerf_dir, num_subfams)
    
    # w0.5 hmms
    w05_dir = kerf_w05_hmm_dir(seed_id, percent_id)
    if os.path.exists(w05_dir):
      os.chdir(w05_dir)
      num_subfams = len(glob.glob(shmm_spec))
      if num_subfams > 0:  hmms_found = True
      if verbose: print "    %d w0.5 shmms" % (num_subfams)
      os.chdir("..")

    # info-share hmms
    infoshare_dir = kerf_infoshare_hmm_dir(seed_id, percent_id)
    if os.path.exists(infoshare_dir):
      os.chdir(infoshare_dir)
      num_subfams = len(glob.glob(shmm_spec))
      if num_subfams > 0:  hmms_found = True
      if verbose: print "    %d info-share shmms" % (num_subfams)
      os.chdir("..")

    os.chdir("..")

  os.chdir(seed_dir)

  # Report sciphy subfamilies
  sciphy_dir = sciphy_hmm_dir(seed_id)
  if os.path.exists(sciphy_dir):
    os.chdir(sciphy_dir)
    num_subfams = len(glob.glob(shmm_spec))
    if num_subfams > 0:  hmms_found = True
    if verbose: print "%d SCI-PHY shmms" % (num_subfams)

  os.chdir(initial_wd)
  return hmms_found

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <seed_id>" % sys.argv[0]
    sys.exit(0)

  seed_id = sys.argv[1]
  report_hmms_for_seed(seed_id, verbose=True)



if __name__ == "__main__":
  main()
