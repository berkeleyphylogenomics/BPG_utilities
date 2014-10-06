#!/usr/bin/env python

from Bio import AlignIO
import sys

h = open(sys.argv[1])

for rec in AlignIO.read(h, 'fasta'):
    rec.seq._data = rec.seq._data.replace('.','-').upper()
    print rec.format('fasta'),
