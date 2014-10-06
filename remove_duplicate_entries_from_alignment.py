#!/usr/bin/env python

import os, sys
from Bio import SeqIO,AlignIO
from Bio.SeqUtils import CheckSum

def remove_duplicate_entries(alignmentfilename):
    # read input alignment
    handle = open(alignmentfilename, 'r')
    alignmentrecords = list(SeqIO.parse(handle, 'fasta'))
    handle.close()

    # get base name based on name of alignment file
    basename = os.path.splitext(alignmentfilename)[0]

    # check alignment for duplicate entries
    num_duplicates = 0
    sequences_of_seguid_of_id = {}
    for record in alignmentrecords:
      id = record.id
      description = record.description
      seq = record.seq.tostring()
      seguid = CheckSum.seguid(seq)
      if id not in sequences_of_seguid_of_id:
        sequences_of_seguid_of_id[id] = {}
      if seguid not in sequences_of_seguid_of_id[id]:
        sequences_of_seguid_of_id[id][seguid] = {}
      if seq in sequences_of_seguid_of_id[id][seguid]:
        num_duplicates += 1
      else:
        sequences_of_seguid_of_id[id][seguid][seq] = description
    # de-dup the input alignment
    if num_duplicates > 0:
      print "Found %d duplicates in %s" % (num_duplicates, alignmentfilename)
      oldalignmentfilename = basename + '_with_dups.afa'
      print "Renaming %s to %s" % (alignmentfilename, oldalignmentfilename)
      os.system("mv %s %s" % (alignmentfilename, oldalignmentfilename))
      print "Writing de-dupped alignment to %s" % alignmentfilename
      f = open(alignmentfilename, "w")
      for id in sequences_of_seguid_of_id:
        for seguid in sequences_of_seguid_of_id[id]:
          for seq in sequences_of_seguid_of_id[id][seguid]:
            f.write(">%s\n" % sequences_of_seguid_of_id[id][seguid][seq])
            f.write("%s\n" % seq)
      f.close()

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <alignment_file_name>" % sys.argv[0]
    sys.exit(0)

  remove_duplicate_entries(sys.argv[1])

if __name__ == '__main__':
  main()
