#!/usr/bin/python
"""
filter_out_dna.py <fasta file with mixed dna and AA sequences>.

- Filters DNA seq out of protein seq file.
- It also filters out sequences exceeding a minimum length
- Filters out seqs with ONLY Xs in it.
- Filters out seqs with ONLY Ns in it.
- Filters out seqs which are mostly DNA (75%).
"""
import os, sys, string, re
from optparse import OptionParser

from BPG_common.sequence import *
from Bio import SeqIO

def main() :
    file = None
    if len(sys.argv) < 2 :
        print 'Need fasta file to process!'
        sys.exit(0)
    else:
        file_path = sys.argv[1]
        if not os.path.exists(file_path):
            print 'File [%s] does not exist!'%file_path
        else:
            file = open(file_path,'r')
    xPattern = re.compile('^[xX]+$')
    nPattern = re.compile('^[nN]+$')
    dna_ratio = 0.75
    #min_seq_length = 40
    min_seq_length = 0
    if file != None:
        for seq_record in SeqIO.parse(file, "fasta") :
            id  = seq_record.id
            name = seq_record.name
            description = seq_record.description
            seq = seq_record.seq.tostring()
#            print '>id: '+id+' name: '+name+' desc: '+description+'\n'
            if not is_dna_or_rna(seq) and len(seq) >= min_seq_length and \
               not xPattern.match(seq) and not nPattern.match(seq) and \
               not is_mostly_dna(seq, dna_ratio):
                print '>' + description + '\n' + seq + '\n'
               

# ------------------------------------------------------------------------------
# Execution of main.
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
# ------------------------------------------------------------------------------

