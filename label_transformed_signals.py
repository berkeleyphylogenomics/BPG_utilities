#!/usr/bin/env python

import sys
from matchmaker.shmm_shmm_lib import *
from numpy import median

def main():
  if len(sys.argv) < 2:
    print "Usage: %s pair" % sys.argv[0]
    sys.exit(0)

  pair = sys.argv[1]
  seedX, seedY = pair.split('_')
  transformed_signals_file = os.path.join(align_dir(seedX, seedY),
                                      "%s_transformed_signals.csv" % pair)
  nr_alignments = set([os.path.split(path)[1] for path in \
              matchmaker_seed_alignment_filenames(seedX, seedY, use_nr=True)])
  f = open(transformed_signals_file)
  header_line = f.readline()
  lines = [line for line in f.readlines() \
            if line.split(',')[0] in nr_alignments]
  f.close()
  q_combined_scores = [float(line.split(',')[1]) for line in lines]
  median_q_combined = median(q_combined_scores)
  sys.stdout.write("Class,")
  sys.stdout.write(header_line)
  for i in range(1,len(lines)):
    if q_combined_scores[i] >= -0.05:
      sys.stdout.write("+,")
      sys.stdout.write(lines[i])
    elif q_combined_scores[i] <= median_q_combined:
      sys.stdout.write("-,")
      sys.stdout.write(lines[i])

if __name__ == '__main__':
  main()
