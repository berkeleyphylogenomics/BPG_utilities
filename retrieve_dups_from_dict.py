#!/usr/bin/python

import os, sys, cPickle

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <runname> <identifier>" % sys.argv[0]
    print "Using the output of "
    print "  make_nr_at_100_with_dict.py <runname> <fasta_file>"
    print "prints all identifiers from <fasta_file> which have the same " \
          + "sequence as <identifier>"""
    sys.exit(0)
  runname = sys.argv[1]
  id = sys.argv[2]
  pklfp = open("%s_dict.pkl" % runname)
  (dups_of_seq, nr_seqs) = cPickle.load(pklfp)
  pklfp.close()
  ids = dups_of_seq[id]
  for id in ids:
    print id
        
if __name__ == '__main__':
  main()

