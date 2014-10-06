#!/bin/env python

import sys
from matchmaker.shmm_shmm_lib import *

pair = sys.argv[1]
(seedX, seedY) = pair.split("_")

print "Retrieving list of all alignment files"
sys.stdout.flush()
afiles = matchmaker_seed_alignment_filenames(seedX, seedY)
n = len(afiles)
print "%d alignment files found" % n
sys.stdout.flush()

i=0
for filename in afiles:
  compress_alignment_file(filename)
  i = i+1
  if i%10 == 0: sys.stdout.write(".")
  sys.stdout.flush()
print "Done."
