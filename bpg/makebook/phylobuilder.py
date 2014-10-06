#!/usr/bin/python

import os
import sys
import string
import re
import time
from optparse import OptionParser
from Bio import SeqIO

"""
phylobuilder.py -> implements the phylobuilder pipeline, which is a way to build
protein family 'books' for the phylofacts resource.
"""


"""
this function just formats a time interval into minutes and seconds,
so we can display the elapsed time for users
"""
def getTimeStr(starttime, endtime):
    elapsed_secs = endtime - starttime

    if elapsed_secs > 60:
        mins = int(elapsed_secs / 60)
        elapsed_secs = elapsed_secs - (mins * 60)

        if mins == 1:
            return '[1 min, %.2f secs]' % elapsed_secs

        return '[%d mins, %.2f secs]' % (mins, elapsed_secs)

    return '[%.2f secs]' % elapsed_secs
# end getTimeStr



"""
Prints a simple 'build-notes' file, which is needed to supply some of the information
for the 'family' table in the database. The build-notes file basically contains
information about how the family was built, which parameter values were used, etc.

This is stored in the file: phylobuilder.buildnotes
"""
def printBuildNotes(ghg):
    handle = open('phylobuilder.buildnotes', 'w')

    if ghg:
        build_method = 'global-global'

    else:
        build_method = 'global-local'

    print >>handle, 'gathering_method FP'
    print >>handle, 'build_clustering_method FLOCK'
    print >>handle, 'build_database_source UniProt'
    print >>handle, 'build_alignment_notes Removed columns with >70% gaps'
    print >>handle, 'family_specific_evalue_criterion 1.0e-5'
    print >>handle, 'family_specific_sw_method %s' % build_method
    print >>handle, 'status normal'
    print >>handle, 'author bpg'
    print >>handle, 'notes normal library build process'

    handle.close()
# end printBuildNotes



"""
gets homologs using the FlowerPower program. input can be either a single protein
sequence 'seed' or a seed alignment. The homologs can be either global-global or
global-local (global to seed, local to gathered homologs).

Resulting alignment should be in final.a2m
"""
def runFlowerPower(inputfile, fromseed, ghg):    
    # run flowerpower from a seed sequence
    if fromseed:
        myfp = 'flowerpower.pl -i %s' % inputfile

    # run flowerpower from an alignment
    else:
        
        # build universe first
        os.system('getUniverseFromMSA.py %s > getUniverseFromMSA.log' % inputfile)
        myfp = 'flowerpower.pl -a %s --msa_profile 1 -u universe.fa' % inputfile

    myfp += ' -n 10 --tempcheck 1'

    if ghg:
        myfp += ' --mode global'

    else:
        myfp += ' --mode glocal'    

    myfp += ' &> flowerpower.log'
    os.system(myfp)

    os.system('perl -pi -e "s/\'//g" final.a2m')
# end runFlowerPower



"""
Aligns sequences using mafft. Uses the number of sequences in the protein family to
determine which algorithm to run. Uses most accurate algorithm if <= 200 sequences -
as per authors' recommendations. If <= 1000 sequences, uses next most accurate
algorithm. Otherwise, lets mafft decide which algorithm to run.

Alignment ends up in phylobuilder.finalaln
"""
def alignSequences(ntax):
    if ntax <= 200:
        mafftcmd = 'mafft --maxiterate 1000 --localpair'

    elif ntax <= 1000:
        mafftcmd = 'mafft --maxiterate 2'

    else:
        mafftcmd = 'mafft'

    mafftcmd += ' phylobuilder.seqs > phylobuilder.finalaln 2> mafft.log'

    os.system(mafftcmd)
# end alignSequences



"""
Runs the build-family pipeline, which builds the phylogenetic tree, pfam, pdb,
etc. Makes a quick decision about which tree-building method to use. If fewer
than 500 sequences, uses ML,NJ. IF fewer than 2000, uses NJ with bootstrap. Otherwise,
uses just NJ without bootstrapping.

Returns a tuple of (NJ,ML), which is true iff that method was run.
"""
def runBuildFamily(ntax, forceml):
    donj = True
    doml = False

    if forceml or ntax <= 500:
        buildfam = 'buildFamily.py --nj --njboot --ml phylobuilder.finalaln'
        doml = True

    elif ntax <= 2000:
        buildfam = 'buildFamily.py --nj --njboot --noml phylobuilder.finalaln'

    else:
        buildfam = 'buildFamily.py --nj --nonjboot --noml phylobuilder.finalaln'

    os.system(buildfam)

    if ntax < 4:
        donj = False
        doml = False

    return (donj, doml)
# end buildFamily



"""
main routine. reads user input, then

1. prints build-notes file
2. runs flowerpower
3. re-aligns sequences using mafft (if global-global option is set)
4. runs build-family (tree, pdb, pfam, etc)
5. runs input-family to put information in the db
6. runs find-orthologs to input ortholog-id stuff into db
"""
def main():
    begintime = time.time()

    myusage = """phylobuilder.py [options] -i seed.fa
  (to build a family from a single seed sequence)
 OR
phylobuilder.py [options] -a alignment.fa
  (to build a family from a seed alignment)
"""
    parser = OptionParser(usage=myusage)

    parser.add_option('--global', action='store_true', dest='global', default=True, help='build family using global homology (default)')
    parser.add_option('--local', action='store_true', dest='local', default=False, help='build family using local homology (global-local)')
    parser.add_option('-i', dest='seedfile', default='NONE', help='input seed file')
    parser.add_option('-a', dest='alnfile', default='NONE', help='input alignment file')
    parser.add_option('--ml', action='store_true', dest='forceml', default=False, help='force construction of ML tree (default: ML tree built only if family has 500 sequences or fewer)')

    (options, args) = parser.parse_args()

    ghg = True

    if options.local:
        ghg = False

    if options.seedfile == 'NONE' and options.alnfile == 'NONE':
        parser.error('either seedfile or alignmentfile required')

    fromseed = True

    if options.seedfile != 'NONE':
        inputfile = options.seedfile

    elif options.alnfile != 'NONE':
        fromseed = False
        inputfile = options.alnfile

    forceml = options.forceml

    ### end user input parsing ###

    # scrape "'" from any input file, as this will muck up the DB entries
    os.system('perl -pi -e "s/\'//g" %s' % inputfile)

    # print build notes file

    printBuildNotes(ghg)

    # get initial alignment using flowerpower

    starttime = time.time()

    print 'running FlowerPower to get initial aligned protein family...',
    sys.stdout.flush()

    runFlowerPower(inputfile, fromseed, ghg)

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)
    sys.stdout.flush()

    # do secondary alignment of sequences (only for global homology)

    ntax = 0

    if ghg:
        starttime = time.time()
        print 'aligning sequences...',
        sys.stdout.flush()

        alnoutf = open('phylobuilder.seqs', 'w')
        handle = open('final.a2m', 'r')

        for seq in SeqIO.parse(handle, 'fasta'):
            myseq = seq.seq.tostring()
            myseq = myseq.replace('.', '')
            myseq = myseq.replace('-', '')
            myseq = string.upper(myseq)

            print >>alnoutf, '>%s' % seq.description
            print >>alnoutf, myseq

            ntax += 1

        handle.close()
        alnoutf.close()

        alignSequences(ntax)

        endtime = time.time()
        print 'done. %s' % getTimeStr(starttime, endtime)
        sys.stdout.flush()   
 
    else:
        # need to do secondary alignment for local books too, but
        # this is a bit trickier. We need to first find the 
        # 'aligned region,' deleting the (potentially large) N- and
        # C-terminal unaligned portions of each sequence, and
        # recording the 'RANGE=[xx..yy]' 
        starttime = time.time()
        print 'aligning matching sub-sequences...',
        sys.stdout.flush()

        alnoutf = open('phylobuilder.seqs', 'w')
        handle = open('final.a2m', 'r')

        mypatt = re.compile('^([a-z.]*)(.+?)[a-z.]*$')

        for seq in SeqIO.parse(handle, 'fasta'):
            myseq = seq.seq.tostring()
            
            # trim N- and C-terminal unaligned regions
            match = mypatt.match(myseq)

            Nunaligned = match.group(1)
            Nunaligned = Nunaligned.replace('.', '')
            mybeginindex = len(Nunaligned)

            myseq = match.group(2)

            myseq = myseq.replace('.', '')
            myseq = myseq.replace('-', '')
            myseq = string.upper(myseq)

            myendindex = mybeginindex + len(myseq)

            # print unalinged sequence
            print >>alnoutf, '>%s' % seq.description + ' RANGE=[%d..%d]' % (mybeginindex, myendindex)
            print >>alnoutf, myseq

            ntax += 1

        handle.close()
        alnoutf.close()

        # align alignable regions
        alignSequences(ntax)

        endtime = time.time()
        print 'done. %s' % getTimeStr(starttime, endtime)
        sys.stdout.flush()

    # make some decisions about if we want an ML tree, bootstrap NJ, and/or
    # Bayesian tree, then run buildFamily

    (didnj, didml) = runBuildFamily(ntax, forceml)

    # run inputFamily    

    cmd = 'inputFamily.py seed.fa phylobuilder.finalaln'
    os.system(cmd)

    # get family "book" id from inputFamily
    handle = open('phylobuilder.familyid', 'r')
    bookid = handle.readline().strip()
    handle.close()

    # find orthologs

    print 'identifying orthlologs...',
    sys.stdout.flush()

    time1 = time.time()

    if didnj:
        cmd = 'find_orthologs --book %s --tree nj > find_orthologs.nj.log' % bookid
        os.system(cmd)

    if didml:
        cmd = 'find_orthologs --book %s --tree ml > find_orthologs.ml.log' % bookid
        os.system(cmd)
    
    cmd = 'find_orthologs --book %s --tree sciphy > find_orthologs.sciphy.log' % bookid
    os.system(cmd)

    time2 = time.time()

    print 'done. %s' % getTimeStr(time1, time2)

    finaltime = time.time()
    print 'phylobuilder finished. %s' % getTimeStr(begintime, finaltime)
    sys.stdout.flush()
# end main



if __name__ == '__main__':
    main()

