#!/usr/bin/python 
import os, glob, sys, re
from matchmaker.shmm_shmm_lib import *
from matchmaker.build_model import *

"""
For each of the top 200 alignments, collect the TSVMod score.
"""
# takes about 8 seconds per model on ohana master node
def collect_TSVMod_score(alignment_filename, pair, score_csv_file):

    if alignment_filename.startswith("/"):
      # we have full pathname for nonmatchmaker alignment
      (seedX, seedY) = pair.split("_")
      fulldirname = os.path.dirname(alignment_filename)
      alignment_filename = os.path.basename(alignment_filename)

    else:
      # otherwise, this is a matchmaker alignment filename
      (seedX, hmmX, seedY, hmmY, filename_tail) = \
	 split_alignment_file_basename(alignment_filename)
      fulldirname = alignment_dir(alignment_filename)

    initial_dir = os.getcwd()
    os.chdir(fulldirname)

    # write scores to master list
    score_filename = "%s.pred_rmsd" % seedY
    if os.path.exists(score_filename):
      score_file = open(score_filename, "r")
      pred_rmsd = float(score_file.readline().strip())
      score_file.close()
    else:
      pred_rmsd = 9999.0

    score_filename = "%s.pred_no35" % seedY
    if os.path.exists(score_filename):
      score_file = open(score_filename, "r")
      pred_no35 = float(score_file.readline().strip())
      score_file.close()
    else:
      pred_no35 = 9999.0

    score_filename = "model_quality_results.csv"
    if os.path.exists(score_filename):
      score_file = open(score_filename, "r")
      score_file.readline()
      try:
        dope = float(score_file.readline().split(",")[1])
      except:
        dope = 9999.0
      score_file.close()
    else:
      dope = 9999.0


    score_csv_file.write("%s, %f, %f, %f\n" % (alignment_filename,
            pred_rmsd, pred_no35, dope))
    
    os.chdir(initial_dir)

def main():

  if len(sys.argv) < 1:
    print "Usage: %s <pair> [<max_models>]" % sys.argv[0]
    sys.exit(0)
  pair = sys.argv[1]
  (seedX, seedY) = pair.split("_")
  if len(sys.argv) > 2: max_models = int(argv[2])
  else: max_models = 200
  score_csv_filename = os.path.join(align_dir(seedX, seedY), "tsvmod_rmsd.csv")
  score_csv_file = open(score_csv_filename, "w")

  # Use list of top alignments from sec-struct/YL ranking
  for alignment_filename in  \
       get_top_ranked_and_special_alignment_filenames(pair, max_models):
      collect_TSVMod_score(alignment_filename, pair, score_csv_file)
  score_csv_file.close()

  # sort results
  cmd = 'sort  -n -t "," --key=2 %s > %s.sorted' % (score_csv_filename,
      score_csv_filename)
  os.system(cmd)

if __name__ == "__main__":
  main()
