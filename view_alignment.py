#!/bin/env python
"""
Input: an alignment filename, either complete path or basename
Output: alignment written to stdout
"""

import sys, os
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print >> sys.stderr, "Usage: %s <alignment_path_or_basename>" % sys.argv[0]
  sys.exit(1)

def display_alignment(filename):
  (seed1, seq1), (seed2, seq2) = read_alignment_file(filename)
  print ">" +  seed1
  print seq1
  print ">" +  seed2
  print seq2

filepath = sys.argv[1]
filename = os.path.basename(filepath)
if filename == filepath:           
  if os.path.exists(filename):
    display_alignment(filename)
  else:
    display_alignment(os.path.join(alignment_dir(filename), filename))
else:
  display_alignment(filepath)
