#!/usr/bin/env python
'''Trim a2m formatted data'''
from Bio import SeqIO, Seq
import re

_trim_re = re.compile(r'^[a-z]*([A-Za-z-]*?[A-Z-])[a-z]*$')

def trim_a2m_record(rec, indices=False):
    '''Trim an a2m SeqRecord object.

    Keyword arguments:
    rec -- Seqrecord object
    indices -- if True, add indices (default False)

    '''

    # replace/remove characters and get start and stop indices
    seq_match = _trim_re.match(rec.seq.tostring().replace('.',''))
    rec.seq = Seq.Seq(seq_match.group(1).replace('-','').upper())

    # get the fasta string and optionally add indices
    if indices is True:
        start, rough_stop = seq_match.span(1)
        stop = rough_stop - seq_match.group(1).count('-')
        rec.id = rec.id.replace(rec.id,
            '%s_%s_%s' % (rec.id, start+1, stop), 1)
    
def _trim_a2m_gen(infile, indices=True):
    for rec in SeqIO.parse(infile, 'fasta'):
        trim_a2m_record(rec, indices)
        yield rec

def trim_a2m(infile, outfile=None, indices=False):
    '''Trim a2m records from an input file handle.

    Keyword arguments:
    infile -- input a2m file handle
    outfile -- output file handle; if None, return a SeqRecord generator
    indices -- if True, add indices (default: False)

    '''

    # if no outfile, return a generator
    if outfile is None:
        return _trim_a2m_gen(infile, indices)

    # write FASTA strings to outfile
    for rec in SeqIO.parse(infile, 'fasta'):
        trim_a2m_record(rec, indices)
        outfile.write(rec.format('fasta'))


if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile',
        help='input A2M file (default STDIN)')
    parser.add_option('-o', '--outfile', dest='outfile',
        help='output FASTA file (default STDOUT)')
    parser.add_option('-x', '--indices', dest='indices',
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

    trim_a2m(infile, outfile, options.indices)

    infile.close()
    outfile.close()
