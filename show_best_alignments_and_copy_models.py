#!/usr/bin/python 
import os, glob, sys, re
from matchmaker.shmm_shmm_lib import *

"""
Show the top N alignments for a pair,
by best TSVmod pred_rmsd score.
Copy the corresponding models into the $CASP_DIR subdirectory
for the target.

Shell variable $CASP_DIR must be defined.
"""

def main():

  if len(sys.argv) < 1:
    print "Usage: %s <pair> [<max_alignments>]" % sys.argv[0]
    sys.exit(0)
  pair = sys.argv[1]
  (seedX, seedY) = pair.split("_")
  if len(sys.argv) > 2: max_alignments = int(sys.argv[2])
  else: max_alignments = 5

  (alignment_files, model_files, alignment_csv_lines) = \
      top_models_by_pred_rmsd(pair, max_alignments)

  target_name = "T0" + seedY[1:]
  target_dir = os.path.join("$CASP_DIR", target_name)

  num_alignments = len(alignment_files)
  for i in range(0, num_alignments):
    alignment_file_basename = os.path.basename(alignment_files[i])
    alignment_file_basename_no_ext = os.path.splitext(alignment_file_basename)[0]
    # print alignment name and scores
    print alignment_csv_lines[i].strip()
    # print alignment
    if os.path.exists(alignment_files[i]):
      """
      7/16/08: replace this code with fn that reads compressed alignment files
      afile = open(alignment_files[i], "r")
      lines = afile.readlines()
      print lines[1].strip()
      print lines[3].strip()
      """
      (seed1, seq1), (seed2, seq2) = read_alignment_file(alignment_files[i])
      print seq1
      print seq2
      print ""
    else:
      print "%s does not exist." % alignment_files[i]
    # copy model into target directory with descriptive name
    new_model_filename = alignment_file_basename_no_ext + ".pdb"
    new_model_pathname = os.path.join(target_dir, new_model_filename)
    if os.path.exists(model_files[i]):
      cmd = "cp %s %s" % (model_files[i], new_model_pathname)
      os.system(cmd)
    else:
      print "%s does not exist." % model_files[i]

if __name__ == "__main__":
  main()
