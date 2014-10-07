""" This is a file for dealing with our alignment software.  Similar to hmm.py for
our hmm functions. """

import subprocess
from pfacts003.fatcat.consts import *
from Bio import AlignIO
from cStringIO import StringIO

def mafft(fasta, options=[]):
    ''' wrapper for using mafft alignment '''
    args = ["/clusterfs/ohana/software/bin/mafft", "--thread", str(THREADS)] + options + ["/dev/stdin"]
    #print args

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(fasta)[0]

def mask_gappy_columns(msa, maximum_fraction_gaps = 0.3):
    ret_msa = []
    msa = AlignIO.read(StringIO(msa), 'fasta')
    num_rows = len(msa)
    
    for column in range(len(str(msa[0].seq))):
        gappiness = 1.0*msa[:,column].count('-')/num_rows
        # first see if the column is sufficiently ungappy
        if (gappiness < maximum_fraction_gaps):
            if ret_msa:
                ret_msa = ret_msa + msa[:,column:(column+1)]
            else:
                ret_msa = msa[:,column:(column+1)]
    if ret_msa:
        return ret_msa.format('fasta')
    else:
        return ""
