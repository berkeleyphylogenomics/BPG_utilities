#!/usr/bin/python

import os, sys, glob
from Bio import SeqIO
from Bio.SeqUtils import CheckSum
import cPickle
from matchmaker.shmm_shmm_lib import *

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <seedX_id> <seedY_id>" % sys.argv[0]
    sys.exit(0)
  seedX_id = sys.argv[1]
  seedY_id = sys.argv[2]
  all_alignment_filenames = \
     matchmaker_seed_alignment_filenames(seedX_id, seedY_id)
  output_handle = open(alignment_nr_csv_file(seedX_id, seedY_id), "w")
  output_handle.write("Alignment\n")
  nr_alignments={}
  dups_of_alignment = {}
  for alignment_filename in all_alignment_filenames:
    """
    #7/15: replaced below with functions that can read compressed alignment files
    f = open(alignment_filename)
    lines = f.readlines()
    if len(lines) < 4:
      continue
    ali = '*'.join([lines[1],lines[3]])
    """
    (seed1, seq1), (seed2, seq2) = read_alignment_file(alignment_filename)
    ali = '*'.join([seq1, seq2])
    the_seguid = CheckSum.seguid(ali)
    if the_seguid in nr_alignments and nr_alignments[the_seguid] != None:
      if ali in nr_alignments[the_seguid] \
          and nr_alignments[the_seguid][ali] != None:
        nr_alignments[the_seguid][ali].add(alignment_filename)
        dups_of_alignment[alignment_filename] = nr_alignments[the_seguid][ali]
      else:
        nr_alignments[the_seguid][ali] = set([alignment_filename])
        output_handle.write("%s\n" % alignment_filename)
        dups_of_alignment[alignment_filename] = nr_alignments[the_seguid][ali]
    else:
      nr_alignments[the_seguid] = { ali: set([alignment_filename]) }
      output_handle.write("%s\n" % alignment_filename)
      dups_of_alignment[alignment_filename] = nr_alignments[the_seguid][ali]

  pklfp = open(os.path.join(align_dir(seedX_id, seedY_id), 
                  "%s_%s_alignment_dict.pkl" % (seedX_id, seedY_id)), "w")
  cPickle.dump((dups_of_alignment, nr_alignments), pklfp)
  pklfp.close()
        
if __name__ == '__main__':
  main()
