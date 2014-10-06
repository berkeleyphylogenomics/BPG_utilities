#!/bin/env python
"""
Input: an alignment filename
Output: directory containing the alignment file
"""

import sys
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print >> sys.stderr, "Usage: %s <alignment_basename>" % sys.argv[0]
  sys.exit(1)

filename = sys.argv[1]
print  alignment_dir(filename)
