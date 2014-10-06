#!/usr/bin/env python

import os, sys
from matchmaker.shmm_shmm_lib import *

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <seed>" % sys.argv[0]
    sys.exit(0)

  seed = sys.argv[1]
  intrepid_dir = os.path.join(cur_dir(), seed, 'INTREPID')
  make_dir_exist(intrepid_dir)
  os.chdir(intrepid_dir)
  os.environ['DATA_DIR'] \
      = '/home/ruchira/src/research/activesite_predictions/intrepid/data'
  f = open('config.txt', 'w')
  f.write("msa_file ../%s.a2m\n" % seed)
  f.write("tree_file ../%s.nj\n" % seed)
  f.write("sequence_id %s\n" % seed)
  f.close()

  status = os.system("intrepid.pl config.txt > intrepid.out")
  if status != 0:
    print "INTREPID exited with nonzero status %d" % status

if __name__ == '__main__':
  main()
