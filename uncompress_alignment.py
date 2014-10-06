#!/bin/env python
"""
Input: an alignment filename, either complete path or basename
Output: file contents, if compressed, replaced with ascii format
"""

import sys, os
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print >> sys.stderr, "Usage: %s <alignment_path_or_basename>" % sys.argv[0]
  sys.exit(1)

filepath = sys.argv[1]
filename = os.path.basename(filepath)
if filename == filepath:           
  if os.path.exists(filename):
    uncompress_alignment_file(filename)
  else:
    uncompress_alignment_file(os.path.join(alignment_dir(filename), filename))
else:
  uncompress_alignment_file(filepath)
