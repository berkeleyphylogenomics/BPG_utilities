#!/usr/bin/python
# Use w0.5 to build HMMs for subfamilies generated
# by kerf (cutting family tree at a certain % ID)

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

# process command line
if len(sys.argv) < 3:
  print "Usage: %s <percent_id> <seed_id>" % sys.argv[0]
  sys.exit(0)
percent_id = int(sys.argv[1])
seed_id = sys.argv[2]

# go to the dir where the kerf subfams are stored
subfamdir = kerf_subfam_dir(seed_id, percent_id)
os.chdir(subfamdir)

# get the list of subfams by listing their fasta seq files
# subfam names are of the form <seed>.kerf<N>.sf<I>
kerf_subfamily_files = glob.glob("%s.kerf%d.sf*.fa" % (seed_id, percent_id))
num_subfams = len(kerf_subfamily_files)
  
# create a dir for the output HMMs, if it doesn't yet exist, and go there
hmmdir = kerf_w05_hmm_dir(seed_id, percent_id)
script_out_filename = os.path.join(hmmdir,"shmm.out")


def w05_already_run(percent_id, seed_id):
  """ See if a .mod file has been created for every .fa file.

      Must be run from dir with hmms.
  """
  hmm_file_spec = "%s.*.sf*.mod" % (seed_id)
  hmm_file_count = len(glob.glob(hmm_file_spec))
  return (hmm_file_count == num_subfams)


make_dir_exist(hmmdir)
os.chdir(hmmdir)
  
if w05_already_run(percent_id, seed_id):
  print "w0.5 shmms already created for kerf%d on %s" \
     % (percent_id, seed_id)
else:

  # remove any .mod files left from an incomplete previous run
  hmm_file_spec = "%s.*.sf*.mod" % (seed_id)
  if len(glob.glob(hmm_file_spec)) > 0:
    cmd = "/bin/rm %s" % (hmm_file_spec)
    os.system(cmd)

  # run w0.5 to build an HMM for every kerf subfamily
  for kerf_subfam_fa in kerf_subfamily_files:
    #take off trailing .fa to get basename for hmm output file
    kerf_subfam = os.path.splitext(kerf_subfam_fa)[0]
    outfile = "%s.mod" % kerf_subfam
    if not (os.path.exists(outfile)):
      # add full pathname prefix to kerf_subfam_fa
      kerf_subfam_fa = os.path.join(subfamdir, kerf_subfam_fa)
      cmd = "w0.5 %s %s >& %s" % (kerf_subfam_fa, outfile,
        script_out_filename)
      os.system(cmd)
  
  # add HMM creation method to alignment filename
  mmv_cmd = 'mmv -d "%s.*.mod" "%s.w05#1.mod"'  % (seed_id, seed_id)
  print mmv_cmd
  os.system(mmv_cmd)
