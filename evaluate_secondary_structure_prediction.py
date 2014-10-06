#!/usr/bin/python

"""
Input:
 -- secondary structure for a protein (given, e.g., by DSSP)
 -- predicted secondary structure for same protein
Both should be in PSIPRED format, in separate files, using the
three-state alphabet.

Output: fraction identity between the two strings.
"""

import sys

def read_ss(ss_file):
  f = open(ss_file)
  lines = f.readlines()
  f.close()
  ss = {}
  for line in lines:
    line = line.rstrip()
    if len(line) > 0 and line[0] != '#':
      fields = line.split()
      ss[fields[0]] = fields[2]
  return ss

if len(sys.argv) < 3:
  print "Usage: " \
    + "%s <dssp_secondary_structure> <predicted_secondary_structure>" \
    % sys.argv[0]
  print """Both the <dssp_secondary_structure> file and the
            <predicted_secondary_structure> file should be in
            PSIPRED VFORMAT."""
  sys.exit(0)
true_ss = read_ss(sys.argv[1])
predicted_ss = read_ss(sys.argv[2])
if len(true_ss.keys()) != len(predicted_ss.keys()):
  print "Error: # of residues differs between the two files"
  print "%d in %s vs %d in %s" \
    % (len(true_ss.keys()), sys.argv[1],
        len(predicted_ss.keys()), sys.argv[2])
  sys.exit(1)
num_residues = float(len(true_ss.keys()))
num_residues_agreeing = 0
for residue_num in true_ss.keys():
  if true_ss[residue_num] == predicted_ss[residue_num]:
    num_residues_agreeing = num_residues_agreeing + 1
q3 = 100 * num_residues_agreeing / num_residues
print "Q3: %0.1f" % q3
