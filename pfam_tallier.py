#!/usr/bin/env python

import sys
import subprocess
from Bio.Blast import NCBIXML
from Bio import SeqIO
import StringIO
import re

proc = subprocess.Popen(['blastall', '-d', 'oldbpg_pfam', '-p', 'blastp', '-i', sys.argv[1], '-m', '7', '-v', '1', '-b', '1', '-e', '0.005'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

f = open(sys.argv[1])
gene_lengths = dict(
    (seq.id, len(str(seq.seq))) for seq in SeqIO.parse(f, 'fasta')
)
f.close()

f = open(f.name.replace('.fasta', '.results'), 'a')

for item in NCBIXML.parse(StringIO.StringIO(proc.communicate()[0])):
    aln = item.alignments[0]
    hsp = aln.hsps[0]
    f.write('\t'.join(map(str, [
        # query id
        item.query,
        # query length
        gene_lengths[item.query],
        # alignment length
        aln.length,
        # query start
        hsp.query_start,
        # query end
        hsp.query_end,
        # subject start
        hsp.sbjct_start,
        # subject end
        hsp.sbjct_end,
        # e value
        hsp.expect,
        # identities/aln_length
        hsp.identities/float(hsp.align_length),
        # hit id
        aln.title,
    ])) + '\n')
    del gene_lengths[item.query]

for key, val in gene_lengths.items():
    f.write('\t'.join(map(str, [key, val])) + '\n')

f.close()
