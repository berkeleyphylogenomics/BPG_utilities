#!/usr/bin/env python

import pickle
import subprocess
import string
import re
import os
import tempfile
from Bio import AlignIO
import glob
import shutil

_aligners = ['mafft', 'muscle']

def satchmojs(infasta, joinlonger=False, dir=None, aligner='mafft',
    iterations=None, kerf_cut=None, minaff=-3.0, out_fasta=None,
    out_lengths_newick=None, out_newick=None, alignment=False, tree=None):

    # check inputs
    if aligner not in _aligners:
        raise Exception("'%s' is not a valid aligner" % aligner)
    elif dir is not None:
        if os.path.exists(dir):
            raise Exception('dir already exists')
        else:
            os.mkdir(dir)
    
    if alignment:
        aligner = 'useralignment'
    i = infasta
    
    tmpdir = dir or tempfile.mkdtemp(
        prefix='satchmojs',
    )
    olddir = os.getcwd()
    os.chdir(tmpdir)

    infile = 'input.fasta'
    aligner_fasta = '%s.fasta' % aligner
    aligner_stockholm = '%s.stockholm' % aligner
    quicktree_newick = 'quicktree.newick'
    satchmo_input = 'satchmo_input'
    subfam_tmp = 'subfam.tmp'
    satchmo_fasta = 'satchmo.fasta'
    satchmo_pretty_fasta = 'satchmo_pretty.fasta'
    satchmo_smo = 'satchmo.smo'
    satchmo_newick = 'satchmo.newick'
    satchmo_phylip = 'satchmo.phylip'
    satchmo_lengths_newick = 'satchmo_lengths.newick'

    
    o = open(infile, 'w')
    # alnstr is used only if the input was an alignment
    # do not replace seq names if a tree was submitted
    # ...there is some ugly control logic going on in here
    #    to make input trees/alignments work
    alnstr = ''
    seq_dict = {}
    counter = 0
    for line in i:
        if line.startswith('>'):
            if tree is None:
                key = 'SEQ%06d' % counter
                seq_dict[key] = line
                keyed = '>%s\n' % key
            else:
                keyed = line
            if alignment:
                alnstr += keyed
            o.write(keyed)
            counter += 1
        else:
            # write the input as unaligned, if it was an alignment
            if alignment:
                alnstr += line
                o.write(line.upper().replace('.','').replace('-',''))
            else:
                o.write(line)
    i.close()
    o.close()

    if tree is None:
        o = open('seq_dict.pkl', 'w')
        pickle.dump(seq_dict, o)
        o.close()


    o = open(aligner_fasta, 'w')
    if not alignment:
        # align
        e = open('%s.err' % aligner, 'w')
        ret = subprocess.call([
            # program
            aligner,
            # iteration argument
            {'mafft': '--maxiterate', 'muscle': '-maxiters'}[aligner],
            # iteration value
            str(iterations is None and \
                {'mafft': 3, 'muscle': 16}[aligner] or iterations),
            # submitted file
            infile,
        ], stdout=o, stderr=e)
        e.close()
        if ret:
            o.close()
            if dir is None:
                shutil.rmtree(tmpdir)
            sys.stderr.write('Alignment step failed.\n')
            sys.exit(1)
    else:
        o.write(alnstr)
    o.close()
    
    # check error output

    e = open('%s.err' % aligner, 'w')
    
    i = open(aligner_fasta, 'r')
    o = open(aligner_stockholm, 'w')

    align = AlignIO.read(i, 'fasta')
    # take the minimum value of the list of pwid's
    pwids = [reduce(lambda m,n: 100.0*sum(m)/sum(n),
        # convert the list of match:align pairs into two lists
        zip(*((
            # do the characters match?
            int(x == y and x != '-'),
            # is the column occupied?
            int(x != '-' or y != '-'),
        # loop over all positions
        ) for x,y in zip(str(a.seq), str(b.seq))))) \
        # loop over all pairs of sequences
        for j,a in enumerate(align) for b in align[j+1:]
    ]
    minpwid = min(pwids)
    maxpwid = max(pwids)
    kerf_cut = kerf_cut is not None and kerf_cut or \
        max((35, minpwid+10))
    AlignIO.write([align], o, 'stockholm')

    i.close()
    o.close()

    if tree is None:
        o = open(quicktree_newick, 'w')
        e = open('quicktree.err', 'w')
        ret = subprocess.call(
            ['quicktree', aligner_stockholm],
            stdout=o, stderr=e
        )
        o.close()
        e.close()

        if ret:
            if dir is None:
                shutil.rmtree(tmpdir)
            sys.stderr.write('Tree construction step failed.\n')
            sys.exit(1)
    else:
        quicktree_newick = 'usertree.newick'
        o = open(quicktree_newick, 'w')
        o.write(tree.read())
        o.close()

    o = open('kerf.out', 'w')
    e = open('kerf.err', 'w')
    ret = subprocess.call(
        ['kerf', '-1', str(kerf_cut), quicktree_newick, aligner_fasta],
        stdout=o, stderr=e,
    )
    o.close()
    e.close()

    if ret:
        if dir is None:
            shutil.rmtree(tmpdir)
        sys.stderr.write('KERF cut step failed\n')
        sys.exit(1)
    
    os.mkdir(satchmo_input)
    summary_partition = open(os.path.join(satchmo_input,
        'summary-partition'), 'w')
    summary_partition.write("Header Line\n")

    for unused, f in sorted((int(item.lstrip('subfam').rstrip('.fa')), item) for item in glob.glob('sub*.fa')):

        subfam_file = os.path.splitext(f)[0]
        summary_partition.write("%s\n" % subfam_file)
        os.mkdir(os.path.join(satchmo_input, subfam_file))

        h = open(subfam_tmp, 'w')
        ret = subprocess.call(['removeGappyColumns', '100', f],
            stdout=h)
        h.close()

        if ret:
            if dir is None:
                shutil.rmtree(tmpdir)
            sys.stderr.write('SATCHMO input preparation step failed.\n')
            sys.exit(1)

        h = open(os.path.join(satchmo_input, subfam_file,
            'acceptedseqs.a2m'), 'w')
        ret = subprocess.call(['removeGappyColumns', '-i', '70', subfam_tmp],
            stdout=h)
        h.close()

        if ret:
            if dir is None:
                shutil.rmtree(tmpdir)
            sys.stderr.write('SATCHMO input preparation step failed.\n')
            sys.exit(1)

        os.unlink(subfam_tmp)
        
    summary_partition.close()
    
    o = open('satchmo.out', 'w')
    e = open('satchmo.err', 'w')
    ret = subprocess.call(['satchmo',
        '-satchmo', infile,
        # conditionally set join method
        ] + (joinlonger and ['-joinmethod', 'longer'] or []) + [
        '-bs', satchmo_input+'/',
        '-out', satchmo_fasta,
        '-sv', satchmo_smo,
        '-phy', satchmo_newick,
        '-minaff', str(minaff),
    ], stderr=e, stdout=o)
    e.close()
    o.close()

    if ret:
        if dir is None:
            shutil.rmtree(tmpdir)
        sys.stderr.write('SATCHMO step failed.\n')
        sys.exit(1)
    
    o = open(satchmo_pretty_fasta, 'w')
    e = open('prettyalign.err', 'w')
    ret = subprocess.call(['prettyalign', satchmo_fasta, '-f', '-m0'],
        stdout=o, stderr=e)
    o.close()
    e.close()

    if ret:
        if dir is None:
            shutil.rmtree(tmpdir)
        sys.stderr.write('prettyalign step failed.\n')
        sys.exit(1)

    i = open(satchmo_pretty_fasta)
    o = open(satchmo_phylip, 'w')
    AlignIO.write(AlignIO.parse(i, 'fasta'), o, 'phylip')
    i.close()
    o.close()

    o = open('raxml.out', 'w')
    e = open('raxml.err', 'w')
    ret = subprocess.call(['raxmlHPC',
        '-f', 'e',
        '-t', satchmo_newick,
        '-m', 'PROTGAMMAJTT',
        '-s', satchmo_phylip,
        '-n', 'post_satchmo_raxml',
        '-w', os.path.abspath('addlengths')],
    stdout=o, stderr=e)
    o.close()
    e.close()

    if ret:
        if dir is None:
            shutil.rmtree(tmpdir)
        sys.stderr.write('raxml step failed.\n')
        sys.exit(1)

    o = open(satchmo_lengths_newick, 'w')
    e = open('iblfut.err', 'w')
    ret = subprocess.call(['incorporate_branch_lengths_from_unrooted_tree.py',
        satchmo_newick,
        'addlengthsRAxML_result.post_satchmo_raxml',
    ], stdout=o, stderr=e)
    o.close()
    e.close()

    if ret:
        if dir is None:
            shutil.rmtree(tmpdir)
        sys.stderr.write('incorporate_branch_lengths_from_unrooted_tree.py step failed.\n')
        sys.exit(1)
    
    newick_trans = string.maketrans(":,)(;][", '_'*7)

    def replace_seqs(filename, handler):

        in_file = open(filename, 'r')
        modified = handler(in_file.readlines())
        in_file.close()

        out_file = open(filename, 'w')
        out_file.writelines(modified)
        out_file.close()

    def replace_fasta(iterable):
        defline_re = re.compile(r'^>(SEQ\d{6})$')
        results = []

        for line in iterable:
            replacement = defline_re.match(line)
            if replacement:
                results.append(seq_dict[replacement.groups()[0]])
            else:
                # If, for whatever reason, no SEQ -- leave original line
                results.append(line)

        return results

    def replace_inline_seqs(iterable):
        seq_re = re.compile(r'(SEQ\d{6})')

        results = []
        for line in iterable:
            replacement = seq_re.split(line)
            for item in replacement:
                if seq_re.match(item):
                    results.append(seq_dict[item].split(' ')[0][1:].\
                                    translate(newick_trans))
                else:
                    results.append(item)

        return results
    if tree is None:
        for filename in [satchmo_pretty_fasta, satchmo_fasta, satchmo_smo]:
            replace_seqs(filename, replace_fasta)

        for filename in [satchmo_lengths_newick, satchmo_newick]:
            replace_seqs(filename, replace_inline_seqs)
    
    if not out_fasta:
        h = open(satchmo_fasta)
        ret = h.read()
        h.close()
    else:
        ret = None

    os.chdir(olddir)
    for outfile,tmp in (
            (out_fasta, satchmo_fasta),
            (out_newick, satchmo_newick),
            (out_lengths_newick, satchmo_lengths_newick),
        ):
        if outfile:
            k = open(os.path.join(tmpdir,tmp))
            outfile.write(k.read())
            k.close()
    if not dir:
        shutil.rmtree(tmpdir)

    return ret

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    
    parser = OptionParser(usage="%prog [options] fasta_file", version="%prog 0.2")
    parser.add_option('-a', '--aligner', dest='aligner', default='mafft',
        choices=_aligners,
        help='choose jump start alignment program (Default: mafft)')
    parser.add_option('-d', '--tmpdir', dest='dir', default=None,
        help='directory to save temporary files')
    parser.add_option('-i', '--iterations', dest='iterations', default=None,
        help='number of iterations for aligner (Default: depends on aligner)')
    parser.add_option('-j', '--joinlonger', action='store_true',
        dest='joinlonger', default=False,
        help='join in SATCHMO by longest hmm')
    parser.add_option('-k', '--kerf_cut', dest='kerf_cut', default=None,
        help='KERF cut % (Default: determined dynamically)')
    parser.add_option('-l', '--lengths_newick', dest='out_lengths_newick',
        default=None, help='write newick tree with lengths')
    parser.add_option('-m', '--minaff', dest='minaff', default=-3.0,
        type='float', help='SATCHMO minimum affinity (Default: -3.0)')
    parser.add_option('-n', '--newick', dest='out_newick',
        default=None, help='write newick tree without lengths')
    parser.add_option('-o', '--out', dest='out_fasta',
        default=None, help='write alignment')
    parser.add_option('-p', '--alignment', dest='alignment',
        default=False, action='store_true',
        help='indicates the input fasta is an alignment')
    parser.add_option('-t', '--tree', dest='tree', default=None,
        help='provide a prepared newick tree; tree MUST have names identical to those found in the input fasta and be compatible with this pipeline (obviates tree construction)')
    options, args = parser.parse_args(sys.argv[1:])
    if len(args) != 1:
        parser.error('incorrect number of arguments')
    elif options.tree is not None and not os.path.exists(options.tree):
        parser.error("'%s' could not be read" % options.tree)
    elif not os.path.exists(args[0]):
        parser.error("'%s' could not be read" % args[0])
    
    for outfile in (options.out_fasta, options.out_lengths_newick,
        options.out_newick):
        if outfile and os.path.exists(outfile):
            parser.error("'%s' already exists" % outfile)
    infasta = open(args[0])
    out_fasta = options.out_fasta is not None and \
        open(options.out_fasta, 'w') or None
    out_lengths_newick = options.out_lengths_newick is not None and \
        open(options.out_lengths_newick, 'w') or None
    out_newick = options.out_newick is not None and \
        open(options.out_newick, 'w') or None
    tree = options.tree is not None and open(options.tree) or None

    out = satchmojs(infasta, joinlonger=options.joinlonger, dir=options.dir,
        aligner=options.aligner, iterations=options.iterations,
        kerf_cut=options.kerf_cut, minaff=options.minaff, out_newick=out_newick,
        out_lengths_newick=out_lengths_newick, out_fasta=out_fasta,
        alignment=options.alignment, tree=tree)
    infasta.close()
    if tree is not None:
        tree.close()

    for outfile in (out_fasta,out_lengths_newick,out_newick):
        if outfile:
            outfile.close()
    if out:
        print out
