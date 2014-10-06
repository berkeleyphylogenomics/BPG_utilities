#!/usr/bin/python

# Reduce size of alignment set by clustering, then selecting
# a representative from each cluster OR computing a consensus
# alignment.

import os, sys
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s  <pair_id>" % sys.argv[0]
  sys.exit(0)

# chdir into the alignment directory.
pair_id = sys.argv[1]
(seedX, seedY) = pair_id.split("_")
out_dir= align_dir(seedX, seedY)
os.chdir(out_dir)
print "Changing to %s" % out_dir

# Do the clustering and find the consensus alignments.
# Assumes that the alignments have been made nonredundant
#  and that the uniquified list of alignment files is in a file
#  <pair>_alignments_nr.csv.
# The full pathnames for the consensus alignments  will be listed in a file,
# consensus_alignments.filenames.
cmd = 'consensus_alignments.py --seedpair %s --nr -v --maxClusters=1000000 &> consensus_alignments.out' % \
 (pair_id)

print "Executing %s" % cmd
os.system(cmd)

alignment_filename = alignment_consensus_csv_file(seedX, seedY)

# Count the final number of alignments
print "Number of alignments after clustering/consensus:"
cmd = "wc %s" % alignment_filename
os.system(cmd)
