#!/usr/bin/python

import sys,os,glob
from matchmaker.shmm_shmm_lib import *

"""
Input: a PDB coordinate file and a chain identifier
Output: same coordinate file, but with all ATOM records
  for other chains removed
"""

if len(sys.argv) < 2:
  print >> sys.stderr, "Usage: %s chain_id (use stdin & stdout)"
  sys.exit(1)

for line in sys.stdin:
  if line[0:4] == "ATOM":
    if line[21] == sys.argv[1]:
      print line.strip()
  else:
    print line.strip()
