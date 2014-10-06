#!/usr/bin/env python
'''Trim and unalign a2m formatted data'''
from Bio import AlignIO, Seq
import re

_trim_re = re.compile(r'^([a-z.]*)([A-Za-z.-]*?[A-Z-])[a-z.]*$')

def trim_a2m_record(rec, indices=False):
    '''Trim an a2m SeqRecord object.

    Keyword arguments:
    rec -- Seqrecord object
    indices -- if True, add indices (default False)

    '''

    # trim
    seq_match = _trim_re.match(rec.seq.tostring())
    rec.seq = Seq.Seq(seq_match.group(2))

    # optionally add indices
    if indices is True:
        rough_start, rough_stop = seq_match.span(2)
        sm1_dots = seq_match.group(1).count('.')
        sm2 = seq_match.group(2)
        rec.id = '%s_%s_%s' % (rec.id,
            # start (adding one for sequence-style indexing)
            rough_start - sm1_dots + 1,
            # stop
            rough_stop - sm1_dots - sm2.count('-') - sm2.count('.'))

def trim_a2m_alignment(alignment, indices=False):
    '''Trim a2m records from an input file handle.

    Keyword arguments:
    alignment -- input a2m Alignment object
    indices -- if True, add indices (default: False)

    '''

    for rec in alignment:
        trim_a2m_record(rec, indices)
    return alignment


if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile',
        help='input A2M file (default STDIN)')
    parser.add_option('-o', '--outfile', dest='outfile',
        help='output FASTA file (default STDOUT)')
    parser.add_option('-r', '--indices', dest='indices',
        help="append indices to record id's (NOTE: varies by header format)",
        action='store_true', default=False)
    options, args = parser.parse_args()

    if len(args):
        import warnings
        warnings.warn(
            'unnecessary arguments: %s' % '; '.join(args),
            UserWarning,
        )

    infile = options.infile is not None and open(options.infile) or sys.stdin
    outfile = options.outfile is not None and open(options.outfile, 'w') or \
        sys.stdout

    for rec in trim_a2m_alignment(
            AlignIO.read(infile, 'fasta'), options.indices):
        outfile.write(rec.format('fasta'))

    infile.close()
    outfile.close()
