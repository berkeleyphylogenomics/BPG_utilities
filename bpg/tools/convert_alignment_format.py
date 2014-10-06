#!/usr/bin/env python

from Bio import AlignIO
import sys

if not len(sys.argv) == 5:
    print "Usage: %s <input file> <input format> <output file> <output format>" % sys.argv[0]
    sys.exit(1)

input_handle = open(sys.argv[1])
alignments = AlignIO.parse(input_handle, sys.argv[2])

output_handle = open(sys.argv[3], "w")
AlignIO.write(alignments, output_handle, sys.argv[4])
 
output_handle.close()
input_handle.close()
