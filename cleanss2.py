#!/bin/env python

#  replace commas with tabs in a CSV file

import sys

if len(sys.argv) < 2:
  print >> sys.stderr, "Usage: %s <file>" % sys.argv[0]
  sys.exit(1)

f = open(sys.argv[1], "r")
lines = f.readlines()
for line in lines[1:]:
  line = line.strip().replace(",", "	")
  print line
