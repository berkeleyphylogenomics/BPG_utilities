#!/usr/bin/python
"""
buildFamily takes care of the non-database processing for building a phylofacts
book. This does NOT include building the alignment, which must be done
beforehand. It also does NOT include putting any of the information in the
database (which is handled by inputFamily) or predicting orthologs (which is
done by find_orthologs, and requires the information to be in the database).

buildFamily takes care of building the phylogenetic tree, searching for PFam
domains, homologous PDB structures, transmembrane helices and signal peptides,
calculating some alignment statistics, and identifying subfamilies using
SCI-PHY.
"""

import os
import re
import sys
import time
from datetime import date
from optparse import OptionParser

import cPickle, cStringIO
from Bio import SeqIO, AlignIO
from Bio.SeqUtils import CheckSum
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML

def getTimeStr(starttime, endtime):

    """Returns time interval formatted in minutes and seconds

    This function formats a time interval into minutes and seconds, so we can
display the elapsed time for users
"""
    elapsed_secs = endtime - starttime

    if elapsed_secs > 60:
        mins = int(elapsed_secs / 60)
        elapsed_secs = elapsed_secs - (mins * 60)

        if mins == 1:
            return '[1 min, %.2f secs]' % elapsed_secs

        return '[%d mins, %.2f secs]' % (mins, elapsed_secs)

    return '[%.2f secs]' % elapsed_secs


def inferNJTree(basename, alignmentdict, bootstrap):

    """Builds the Neighbor Joining tree using Quicktime

    Input alignment gets printed to .aln.nj. Raw output tree is in .nj.tre, and
rooted tree is in .nj.rooted.tre

    If desired, we'll do 200 bootstrap reps to calculate clade support.
However, note that mpreroot (ie, retree wrapper) dumps the support values, so
the rooted tree will not have support values (boo hoo).

    We are now using our own midpoint_reroot.py instead of mpreroot when we do
    bootstrap reps.  This may not dump the support values, but this has not
    been tested yet.


"""

    bootstrapreps = 200

    # print alignment (stockholm format)
    njalnfname = basename + '.aln.nj'
    handle = open(njalnfname, 'w')

    for k in alignmentdict.keys():
        print >>handle, '%s %s' % (k, alignmentdict[k])

    handle.close()

    # infer nj tree

    cmd = 'quicktree -in a -out t '

    if bootstrap:
        cmd += '-boot %d ' % bootstrapreps

    njtreefname = basename + '.nj.tre'
    cmd += '%s > %s' % (njalnfname, njtreefname)

    os.system(cmd)

    # root nj tree at midpoint

    # mpreroot (which just calls retree) dumps support values
    # In that case, midpoint_reroot.py instead, which hopefully won't dump the
    # support values.
    # However, midpoint_reroot.py (written in Python) is slower than mpreroot,
    # so if we have a plain neighbor-joining tree, we run mpreroot.
    njrootedtreefname = basename + '.nj.rooted.tre'
    if bootstrap:
      cmd = 'export PYTHONPATH=/clusterfs/ohana/software/lib/python2.6; ' \
          + 'midpoint_reroot.py %s > %s' % (njtreefname, njrootedtreefname)
    else:
      cmd = 'mpreroot %s > %s' % (njtreefname, njrootedtreefname)
    os.system(cmd)



def inferMLTree(basename, alignmentdict, use_fasttree=True):

    """Build Maximum Likelihood Tree using phyml or FastTree.

    With FastTree we use the default paramters, except we use the -gamma option
    to get the approximate Gamma(20) log-likelihood values, which are more
    comparable across runs.

    We assume the JTT+G4 model input alignment is in .aln.ml, output tree is in
.ml.tre, and rooted tree (rooted at midpoint) is in .ml.rooted.tre.

    We calculate support as the approximate likelihood ratio (aLR) for each
node, which is the likelihood of the best tree divided by the likelihood of the
next-best tree.

    Note that mpreroot (ie, retree) dumps all support values, so the support
measures will NOT be in the rooted tree.

    We are now using our own midpoint_reroot.py instead of mpreroot.  This does
    not dump the likelihood values output by FastTree.  It has not been tested
    with phyml output yet.
"""
    
    # print alignment file in phylip format

    mlalnfname = basename + '.aln.ml'
    handle = open(mlalnfname, 'w')

    ntax = len(alignmentdict.keys())
    nchar = 0

    for k in alignmentdict.keys():
        nchar = len(alignmentdict[k])
        break

    print >>handle, '%d %d' % (ntax, nchar)

    for k in alignmentdict.keys():
        print >>handle, '%s %s' % (k, alignmentdict[k])

    handle.close()

    if use_fasttree:
      # Infer Maximum Likelihood tree with FastTree.  Use default parameters,
      # except use -gamma to compute approximate Gamma(20) log-likelihoods,
      # which are comparable across runs (unlike the CAT log-likelihoods).

      cmd = 'FastTree -gamma %s > %s.fasttree.ml.tre 2> %s.fasttree.log' \
            % (mlalnfname, basename, basename)

      os.system(cmd)

      # Root ML tree at midpoint

      cmd = 'export PYTHONPATH=/clusterfs/ohana/software/lib/python2.6; ' \
          + 'midpoint_reroot.py %s.fasttree.ml.tre > %s.fasttree.ml.rooted.tre' \
            % (basename, basename)

      os.system(cmd)
      return

    #
    # infer best-fit protein model using ProtTest. We can use
    # AIC criterion.
    #

    cmd = 'prottest -i %s -o %s.prottest &> prottest.out' % (mlalnfname, mlalnfname)
    os.system(cmd)

    prottestpatt = re.compile('Best model according to AIC: (.+)')

    handle = open(mlalnfname + '.prottest', 'r')

    modelname = ''
    invar = False
    gamma = False
    freq  = False

    for line in handle:
        match = prottestpatt.match(line)

        if match:
            rawmodel = match.group(1).strip()
    
            # print model to output file
            handle2 = open(basename + '.mlmodel', 'w')
            print >>handle2, rawmodel
            handle2.close()

            modelarr = rawmodel.split('+')

            for x in modelarr:
                if x == 'I':
                    invar = True
                elif x == 'G':
                    gamma = True
                elif x == 'F':
                    freq = True
                else:
                    modelname = x

            break

    handle.close()

    # infer ML tree using best-fit model

    cmd  = '/clusterfs/ohana/software/ProtTest2.0/bin/phyml-prottest-linux64'
    cmd += ' -i %s -d aa -q -b -1 -m %s' % (mlalnfname, modelname)
    
    if freq:
        cmd += ' -f e'
    else:
        cmd += ' -f m'

    if invar:
        cmd += ' -v e'
    else:
        cmd += ' -v 0.0'

    if gamma:
        cmd += ' -c 4 -a e'
    else:
        cmd += ' -c 1'

    cmd += ' -s NNI -o tlr &> %s.phyml.log' % basename

    os.system(cmd)

    # change name of tree file

    mltreefname = basename + '.ml.tre'
    os.system('mv %s_phyml_tree.txt %s' % (mlalnfname, mltreefname))

    # root ml tree at midpoint

    mlrootedtreefname = basename + '.ml.rooted.tre'
    cmd = 'mpreroot %s > %s' % (mltreefname, mlrootedtreefname)
    cmd = 'export PYTHONPATH=/clusterfs/ohana/software/lib/python2.6; ' \
        + 'midpoint_reroot.py %s > %s' \
          % (mltreefname, mlrootedtreefname)
    os.system(cmd)



def createHMM(basename, alignmentfilename):

    """Creates a HMM from a sequence alignment, using SAM w0.5

    Creates a hidden markov model (HMM) from a sequence alignment, using the
SAM w0.5 software. Also converts the model to HMMER format.

    SAM model goes in .mod
    HMMER model goes in .hmm

    consensus sequence generated from HMM goes in .con.fa
"""

    # create initial HMM using extra-special w0.5 method

    modfname = basename + '.mod'
    cmd = 'w0.5 %s %s &> w0.5.out' % (alignmentfilename, modfname)
    os.system(cmd)

    # convert binary model to ASCII
    cmd = 'hmmconvert %s -model_file %s' % (basename, modfname)
    os.system(cmd)

    f = open(alignmentfilename)
    alignment = AlignIO.read(f, 'fasta')
    f.close()
    # Convert alignment to Stockholm format
    # Save in a string for later use, as we will insert a line
    s = cStringIO.StringIO()
    AlignIO.write([alignment], s, "stockholm")
    stockholm_lines = s.getvalue().split('\n')
    s.close()

    # Get the first record from the alignment
    for record in alignment:
      break

    # Write the Stockholm format file
    stockholm_fname = basename + '.stockholm'
    f = open(stockholm_fname, "w")
    # Write everything before the // line
    f.write('\n'.join(stockholm_lines[:-2]))
    f.write('\n')
    # Write the #=RF line indicating which columns are aligned
    f.write("#=GC RF ")
    for ch in record.seq:
      if ch.isupper() or ch == '-':
        f.write("x")
      else:
        f.write(".")
    f.write("\n")
    # Write the // line
    f.write('\n'.join(stockholm_lines[-2:]))
    f.write('\n')
    f.close()

    # create initial HMMER3 hmm
  
    hmmfname = basename + '.hmm'
    cmd = 'hmm3build -n %s --hand %s %s &> hmm3build.out' \
          % (basename, hmmfname, stockholm_fname)
    os.system(cmd)

    # generate consensus sequence from HMM

    consensfname = basename + '.con.fa'
    cmd = 'hmm3emit -c %s > %s' % (hmmfname, consensfname)
    os.system(cmd)



def inferPFam(basename):
    """infers PFam domains from the consensus sequence using hmmpfam.

    PFam results go in .pfam

    Note that this step seems to take a long time; speeding it up could help.
"""

    # uses hmmscan to score the consensus sequence against
    # the PFam HMM database    

    seqfname = basename + '.con.fa'
    verbose_pfamfname = basename + "_vs_Pfam-A.hmmscan.out"
    domain_tablefname = basename + "_vs_Pfam-A.domtbl.out"
    cmd = 'hmmscan --cut_ga ' \
        + '-o %s ' % verbose_pfamfname \
        + '--domtblout %s ' % domain_tablefname \
        + '/clusterfs/ohana/external/pfam/current/Pfam-A.hmm ' \
        + '%s ' % seqfname
    os.system(cmd)


def inferTransmembrane(basename):

    """Infers transmembrane helices and signal peptides from consensus sequence

    Infers transmembrane helices and signal peptides from the family consensus
sequence using the Phobius webserver. Note that we need web access to do this.

    results file is .phobius
"""

    # read consensus sequence from file

    seqfname = basename + '.con.fa'

    handle = open(seqfname, 'r')

    seq = ''

    line = handle.readline()
    line = handle.readline().strip()

    while line:
        seq += line
        line = handle.readline().strip()

    handle.close()

    # use phobius webserver to predice transmembrane helices
    # and signal peptides from the consensus sequences

    phobiusfilename = basename + '.phobius'
    cmd = 'ssh ohana curl -F "protseq=%s" -F "format=nog" http://phobius.binf.ku.dk/cgi-bin/predict.pl > %s 2> phobius.out' % (seq, phobiusfilename)
    os.system(cmd)



def inferPDB(basename):

    """Infers homologous PDB structures from the family HMM.

    Score the PDB sequence database with the family HMM.

    Results go in <basename>_vs_PDB.hmmsearch.out and <basename>_vs_PDB.tbl.out
"""
    
    verbose_pdbhit_fname = basename + "_vs_PDB.hmmsearch.out"
    pdbhit_table_fname = basename + "_vs_PDB.tbl.out"
    hmmfname = basename + '.hmm'
    cmd = 'hmm3search '\
          + '-o %s ' % verbose_pdbhit_fname \
          + '--tblout %s ' % pdbhit_table_fname \
          + '%s ' % hmmfname \
          + '/clusterfs/ohana/external/pdb_rcsb_full'
    os.system(cmd)

def inferSubfamilies(basename, alignmentdict):

    """Infers subfamilies using SCI-PHY

    SCI-PHY uses the alignment and prints a bunch of output files.
"""

    # print alignment with internal ids that sci-phy can use

    alnfname = basename + '.aln.sciphy'

    handle = open(alnfname, 'w')

    for k in alignmentdict.keys():
        print >>handle, '>%s' % k
        print >>handle, alignmentdict[k]

    handle.close()

    # use SCI-PHY to infer subfamilies and generate tons of output files

    # why does SCI-PHY suck so much? because it crashes, very often.
    # apparently, no one knows why. but, you have to just keep trying
    # to run it until it happens to want to work this time.

    cmd = 'SCI-PHY %s -i %s -hmm %s -print_tree &> sciphy.out' % (basename, alnfname, basename + '.mod')

    i = 0    

    while i < 5 and os.system(cmd) != 0:
        i += 1

def computeAlignmentConservation(basename, alignmentfilename):
  outfname = basename + '.alignmentconservation.csv'
  os.system('compute_alignment_conservation %s > %s' % (alignmentfilename,
                                                        outfname))

def getAlignmentStatistics(basename, alignmentfilename, astats):

    """Gather alignment statistics on input alignment

    Uses astats to gather a bunch of alignment statistics on the input alignment.  Also prints NR100 (non-redundant alignment at 100%), which is needed to fill in some of the database entries in the family table.

    Results go in .astats and .nr100.a2m, respectively.
"""
    # get astats statistics
    if astats:
      outfname = basename + '.astats'
      os.system('astats %s > %s' % (alignmentfilename, outfname))

    # get NR100 alignment
    os.system('make_nr_at_100_with_dict.py %s %s' % (basename, 
                                                    alignmentfilename) 
              + '>& make_nr_at_100_with_dict.out')

def buildFamily(alignmentfilename, njtree, njbootstrap, mltree, sciphy, astats,
                  fasttree):
    """main routine for buildFamily, runs the pipeline..."""

    starttime = time.time()
    print 'reading input alignment %s...' % alignmentfilename,
    sys.stdout.flush()
    
    # get base name based on name of alignment file
    # this will be used to create names for tree files, etc
    basename = os.path.splitext(alignmentfilename)[0]

    # read input alignment
    handle = open(alignmentfilename, 'r')
    alignmentrecords = list(SeqIO.parse(handle, 'fasta'))
    handle.close()

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
      # read the input alignment again, so it is de-dupped
      handle = open(alignmentfilename, 'r')
      alignmentrecords = list(SeqIO.parse(handle, 'fasta'))
      handle.close()

    # create alignment for treebuilding
    # also create an ID mapping from internal to fasta identifiers
    idmap = {}
    treealignment = {}
    alignmentdict = {}
    id = 1

    for record in alignmentrecords:
        myid = 'SEQ%d' % id
        id += 1
    
        alignmentdict[myid] = record.seq.tostring()

        idmap[myid] = record.description        

        treealignment[myid] = ''

        for c in record.seq.tostring():
            if c == '-' or c.isupper():
                treealignment[myid] += c
    
    # print ID map file (pickle)
    idmapfname = basename + '.idmap'
    handle = open(idmapfname, 'w')
    cPickle.dump(idmap, handle)
    handle.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # print alignment file for NJ and build NJ tree
    if njtree:
        ntaxa = len(alignmentrecords)

        if ntaxa >= 4:
            
            starttime = time.time()
            print 'inferring tree by neighbor joining...',
            sys.stdout.flush()

            inferNJTree(basename, treealignment, njbootstrap)

            endtime = time.time()
            print 'done. %s' % getTimeStr(starttime, endtime)

        else:
            print 'too few sequences (%d) to build neighbor joining tree, skipping.' % (ntaxa)

    # print alignment for ML tree and build ML tree
    if mltree:
        ntaxa = len(alignmentrecords)

        if ntaxa >= 4:

            starttime = time.time()
            print 'inferring tree by maximum likelihood...',
            sys.stdout.flush()

            inferMLTree(basename, treealignment, fasttree)
    
            endtime = time.time()
            print 'done. %s' % getTimeStr(starttime, endtime)

        else:
            print 'too few sequences (%d) to build maximum likelihood tree, skipping.' % (ntaxa)

    # create family-level hmm

    starttime = time.time()
    print 'creating general hidden Markov model...',
    sys.stdout.flush()

    createHMM(basename, alignmentfilename)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # get PFam domains

    starttime = time.time()
    print 'inferring PFam domains...',
    sys.stdout.flush()

    inferPFam(basename)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # get transmembrane and signal peptide predictions

    starttime = time.time()
    print 'inferring transmembrane domains and signal peptides...',
    sys.stdout.flush()

    inferTransmembrane(basename)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # score PDB

    starttime = time.time()
    print 'retrieving homologous PDB structures...',
    sys.stdout.flush()

    inferPDB(basename)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # run SCI-PHY to infer subfamilies
    if sciphy:

      starttime = time.time()
      print 'inferring subfamilies...',
      sys.stdout.flush()

      inferSubfamilies(basename, alignmentdict)

      endtime = time.time()
      print 'done. %s' % getTimeStr(starttime, endtime)

    # compute alignment conservation
    starttime = time.time()
    print 'calculating alignment conservation...',
    sys.stdout.flush() 

    computeAlignmentConservation(basename, alignmentfilename)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # gotta run astats to get bulks of info about who's long, who's short,
    # and who's dating whom

    starttime = time.time()
    print 'calculating alignment statistics...',
    sys.stdout.flush()

    getAlignmentStatistics(basename, alignmentfilename, astats)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    # record build date

    print 'recording build date...',
    sys.stdout.flush()
    
    datefname = basename + '.build_date'
    handle = open(datefname, 'w')
    print >>handle, date.today()
    handle.close()

    print 'done.'

def main():
    """Main routine, so that buildFamily can be run as a script.

    This just parses user input and then runs buildFamily function.
"""

    parser = OptionParser(usage='buildFamily.py [options] alignmentfile')

    parser.add_option('--sciphy', action='store_true', dest='sciphy',
        default=True,
        help="run SCI-PHY automated subfamily identification")
    parser.add_option('--nosciphy', action='store_false', dest='sciphy',
        default=True,
        help="do not run SCI-PHY automated subfamily identification")
    parser.add_option('--astats', action='store_true', dest='astats',
        default=True,
        help="run astats to compute alignment statistics")
    parser.add_option('--noastats', action='store_false', dest='astats',
        default=True,
        help="do not run astats to compute alignment statistics")
    parser.add_option('--nj', action='store_true', dest='nj',
        default=True,
        help='infer tree using neighbor joining')
    parser.add_option('--njboot', action='store_true', dest='njboot',
        default=True,
        help='calculate bootstrap support for neighbor joining tree')
    parser.add_option('--nonjboot', action='store_false', dest='njboot',
        default=True,
        help='calculate bootstrap support for neighbor joining tree')
    parser.add_option('--nonj', action='store_false', dest='nj',
        default=True,
        help='do not infer tree using neighbor joining')
    parser.add_option('--ml', action='store_true', dest='ml',
        default=True,
        help='infer tree using maximum likelihood')
    parser.add_option('--noml', action='store_false', dest='ml', 
        default=True,
        help='do not infer tree using maximum likelihood')
    parser.add_option('--fasttree', action='store_true', dest='use_fasttree',
        default=True,
        help='use FastTree for inferring maximum likelihood trees')
    parser.add_option('--nofasttree', action='store_false',
        dest='use_fasttree',
        default=True,
        help='use PhyML, not FastTree, for inferring maximum likelihood trees')

    (options, args) = parser.parse_args()

    njtree = options.nj
    njboot = options.njboot
    mltree = options.ml
    sciphy = options.sciphy
    astats = options.astats
    fasttree = options.use_fasttree

    if len(args) < 1:
        parser.error('alignment filename argument required')

    alnfname = args[0]
    buildFamily(alnfname, njtree, njboot, mltree, sciphy, astats, fasttree) 


if __name__ == '__main__':
    main()

