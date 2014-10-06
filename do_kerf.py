#!/usr/bin/python

# Cut tree at given %id using an alignment

import os, sys, time, glob
from matchmaker.shmm_shmm_lib import *

def kerf_already_completed(seed_id):
  """
  Must be run from kerf results directory.
  """
  # if .summary file does not exist, return False
  kerf_summary_filename = "kerf.summary"
  if not os.path.exists(kerf_summary_filename): return False
  #print "%s exists" % kerf_summary_filename
  # read last line of summary file and get subfam #
  # if last line just has dashes, no subfams were created
  kerf_summary_file = open(kerf_summary_filename, "r")
  line = kerf_summary_file.readline()
  prev_line = line
  while (line):
    prev_line = line
    line = kerf_summary_file.readline()
  #print "line = %s" % line
  #print "prev_line = %s" % prev_line
  tokens = prev_line.split()
  #print "tokens: ", tokens
  if len(tokens) == 0:
    #print "No tokens"
    return False
  if (tokens[0].startswith("----")):
    summary_count = 0
  else:
    try:
      summary_count = int(tokens[0]) + 1
    except:
      return False
  # count number of subfamily files
  tree_count = len(glob.glob("%s.kerf*.tre" % seed_id))
  msa_count = len(glob.glob("%s.kerf*.fa" % seed_id))
  #print "summary, tree, msa counts" , summary_count, tree_count, msa_count
  # if counts match, return true
  if (summary_count == tree_count) and (tree_count == msa_count):
    return True
  else:
    return False


def main():
  if len(sys.argv) < 3:
    print "Usage: %s <percent_id> <seed_id>" % sys.argv[0]
    sys.exit(0)

  percent_id = int(sys.argv[1])
  seed_id = sys.argv[2]
  root = "%s/%s" % (cur_dir(), seed_id)
  out_dir = os.path.join(root, "kerf%d" % percent_id)
  make_dir_exist(out_dir)
  os.chdir(out_dir)

  if kerf_already_completed(seed_id):
    print "kerf%d previously completed for %s" % (percent_id, seed_id)
  else:
    cmd = '%s "-1" %d ../%s.nj ../%s.a2m' \
            % (kerf_cmd(), percent_id, seed_id, seed_id)
    print "Running:", cmd
    os.system(cmd)
    cmd = 'mmv -d "subfam*.*" "%s.kerf%d.sf#1.#2"' % (seed_id, percent_id)
    print "Running:", cmd
    os.system(cmd)

if __name__ == "__main__":
  main()
