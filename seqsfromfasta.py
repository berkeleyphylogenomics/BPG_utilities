#!/usr/bin/env python

from Bio import SeqIO

def seqpull(h, *args): #should use 'any' in py > 2.3
    return ''.join([seq.format('fasta') for seq in SeqIO.parse(h,'fasta') \
        if sum([seq.id.count(arg) for arg in args])])

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print "%s: get sequences from a fasta file by substring in defline" \
            % sys.argv[0]
        print "USAGE: %s <multiple fasta file> [keywords]" % sys.argv[0]
    else:
        h = open(sys.argv[1])
        print seqpull(h,*sys.argv[2:])
        h.close()
