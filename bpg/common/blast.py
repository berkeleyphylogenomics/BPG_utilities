#!/usr/bin/python

# Blast interface
#

import cPickle
import glob
import os
import re
import string
import sys
import tempfile

import config
from bpg.common import fasta
from bpg.common.legacy_utils import *

BLASTDB = os.environ.get( "BLASTDB" )

Verbose = 0

# -----------------------------------------------------------------

# Returns a dictionary where the keys are the names and the values are
# (sequence, evalue) pairs

def RunBlast(seeds, dbfile, how, biter, maxseqs, maxseqsper, tmpbase = "", cachefile = "", dbdir = BLASTDB) :
    if cachefile and os.access(cachefile, os.F_OK) :
        fp = open(cachefile, "r")
        results = cPickle.load(fp)
        fp.close()
        print "%d blast results loaded from cache %s" % (len(results), cachefile )
        return results
        

    # These should come from standard environment.

    ###os.environ["BLASTDB"] = dbdir
    ###os.environ["BLASTMAT"] = config.BLASTMAT

    if tmpbase :
        tmp = tmpbase
    else :
        tmp = tempfile.mktemp()
    
    results = {}

    num = 1
    for s in seeds :
        if Verbose : print "BLAST %s" % s
        
        infile = "%s.seed%d.blastin" % (tmp, num)
        outfile = "%s.seed%d.blastout" % (tmp, num)
        num += 1
        
        seq = fasta.OnlyMatches(seeds[s])
        file(infile, "w").write(seq)

        # These are the parameters flowerpower wants.
        if len(seq) < 35 :
            matrix = "PAM30"
            gapcost = 9
            alignment = 50
            e = 1000
        else :
            alignment = 1000
            gapcost = 10
            e = 10
            if len(seq) < 50 :
                matrix = "PAM70"
            elif len(seq) <= 85 :
                matrix = "BLOSUM80"
            else :
                matrix = "BLOSUM62"

        outfile_p = outfile + ".pgp"
        outfile_a = outfile + ".all"

        # -o output file
        # -d database
        # -i queryfile
        # -m 9 output Tabular output with comments
        # -j number of passes
        # -h eval threshhold
        # -e expectation value (?)
        # -v number of sequences to show one-line descriptions for
        # -K number of best hits from a region to keep
        # -G gap open cost
        # -M matrix 
        # -b Number of database sequence to show alignments for 
        cmd_pgp = "%sblastpgp -o %s -d %s -i %s -m 9 -j %s -h 0.001 -e %s -v 1000 -K 1000 -G %s -M %s -b %s" % (config.BLASTDIR, outfile_p, dbfile, infile, biter, e, gapcost, matrix, alignment)
        cmd_all = "%sblastall -p blastp -o %s -d %s -i %s -m 9 -e %s -v 1000 -K 1000 -G %s -M %s -b %s" % (config.BLASTDIR, outfile_a, dbfile, infile, e, gapcost, matrix, alignment)

        rr = []
        if how == "blastpgp" or how == "both" :
            if not os.access(outfile_p, os.F_OK) :
                RunCmd(cmd_pgp, "", Verbose)
            elif Verbose :
                print "not running blast, %s found" % outfile_p
            rr = rr + GetResults(outfile_p, dbfile)
            if not tmpbase : os.unlink(outfile_p)

        if how == "blastall" or how == "both" :
            if not os.access(outfile_a, os.F_OK) :
                RunCmd(cmd_all, "", Verbose)
            elif Verbose :
                print "not running blast, %s found" % outfile_a
            rr = rr + GetResults(outfile_a, dbfile)
            if not tmpbase : os.unlink(outfile_a)

        if maxseqsper :
            rr.sort()
            rr = rr[0 : maxseqsper]
            
        # Note that psiblast includes results multiple times for the different
        # iterations so we want to filter them here.
        for (e, n, s) in rr :
            if not results.has_key(n) or results[n][1] > e :
                results[n] = (s, e)

        if Verbose : print "total results = %d" % len(results)
                
        if not tmpbase :
            os.unlink(infile)

    if maxseqs and len(results) > maxseqs :
        rr = []
        for r in results : rr.append((results[n][1], results[n][0], n))
        rr.sort()
        rr = rr[0 : maxseqs]
        results = {}
        for r in rr : results[r[2]] = [r[1], r[0]]
    
    if cachefile :
        fp = open(cachefile, "w")
        cPickle.dump(results, fp, 1)
        fp.close()
        print "%d blast results written to cache %s" % (len(results), cachefile )
        
    return results
    
# -----------------------------------------------------------------

def LoadAlignment(queryname, blastout) :
    defline = ""
    defline1 = ""
    query = []
    sbjct = []
    ret = []
    for line in file(blastout) :
        if len(line) < 2 : continue
        if line[0] == ">" :
            if query :
                ret.append((defline, query, sbjct))
                query = []
                sbjct = []
            defline1 = line[1:].strip()
        elif defline1 :
            if line.strip().startswith("Length =") :
                defline = defline1
                defline1 = ""
                continue

            elif line[1] != " " :
                line = "|" + line[1:].strip()
            else :
                line = " " + line.strip()
            defline1 = defline1 + line
        elif line.startswith("Query: ") :
            query.append(line.split()[1:3])
        elif line.startswith("Sbjct: ") :
            sbjct.append(line.split()[1:3])
    if query :
        ret.append((defline, query, sbjct))

    qseq = {}
    sseqs = {}
    maxp = 0
    for r in ret :
        defline = r[0]
        sseq = {}
        sseqs[defline] = sseq
        nl = len(r[1])
        if len(r[1]) != len(r[2]) : Error("format hosed")
        for i in range(nl) :
            qn = int(r[1][i][0])
            sn = int(r[2][i][0])
            q = r[1][i][1]
            s = r[2][i][1]
            for j in range(len(q)) :
                qc = q[j]
                sc = s[j]
                if qc != '-' :
                    if qseq.get(qn, qc) != qc : Error("mismatch")
                    qseq[qn] = qc
                    sseq[qn] = sc
                    if qn > maxp : maxp = qn
                    qn = qn + 1
    ret2 = {}
    for i in range(1, maxp + 1) :
        qs = ret2.get(queryname, "")
        if not qseq.has_key(i) : Error("gap found in query position %d" % i)
        ret2[queryname] = qs + qseq[i]
        for s in sseqs :
            ss = ret2.get(s, "")
            ret2[s] = ss + sseqs[s].get(i, "-")
    return ret2

# -----------------------------------------------------------------

def GetResults(blastout, dbfile) :
    results = []
    fp = file(blastout)
    reading = 0
    lookup = []
    evals = {}
    for line in fp :
        if line.startswith("#") : continue

        line1 = line.split()
        if len(line1) < 12 : continue

        name = line1[1]
        eval = line1[10]
        
        comp1 = string.join(name.split('|')[0:2], '|')
        evals[comp1] = eval
        lookup.append(comp1)

    if lookup :
        seq = Fastacmd(lookup, dbfile)
        for s in evals :
            got = ""
            for s1 in seq :
                if s1.find(s) >= 0 :
                    got = seq[s1]
                    break
            if not got : raise "Can't find %s in fasta" % s
            results.append((float(evals[s]), s, got))
        
    if Verbose : print "read %s from %s" % (Plural("result", len(results)), blastout)
            
    return results
    
# -----------------------------------------------------------------

def GetResultsFromHumanReadable(blastout, dbfile) :
    results = []
    fp = file(blastout)
    reading = 0
    lookup = []
    evals = {}
    for line in fp :
        if line.startswith(">") :
            break
        
        elif reading :
            try :
                eval = float(line[76:])
            except :
                continue
            name = line[0:64]
            comp1 = string.join(name.split('|')[0:2], '|')
            lookup.append(comp1)
            evals[comp1] = eval
            
        elif line.startswith("Sequences producing significant alignments") :
            reading = 1

    if lookup :
        seq = Fastacmd(lookup, dbfile)
        for s in evals :
            if not seq.has_key(s) : raise "Can't find %s in fasta" % s
            results.append((s, seqs[s], evals[s]))
        
    if Verbose : print "read %d results from %s" % (len(results), blastout)
            
    return results
    
# -----------------------------------------------------------------

def Fastacmd(names, dbfile = None) :
    prog = config.BLASTDIR + "fastacmd"
    tmp = tempfile.mktemp()
    infile = tmp + ".in"
    outfile = tmp + ".fa"
    fp = file(infile, "w")
    for name in names : fp.write(name + '\n')
    fp.close()
    if dbfile :
        cmd = "%s -i %s -d %s" % (prog, infile, dbfile)
    else :
        cmd = "%s -i %s" % (prog, infile)
    if not RunCmd(cmd, outfile, "", -1) :
        seqs = fasta.ReadSequences(outfile)
    else :
        seqs = {}
        if len(names) == 1 :
            Print(cmd)
            Print("Warning: fastacmd can't look up: %s" % names)
            Print(file(outfile).read())
        else :
            # At least one of those GI's didn't work.  Let's
            # look them all up separately and report the ones that failed.
            for name in names :
                seqs2 = Fastacmd([name], dbfile)
                seqs.update(seqs2)
    os.remove(infile)
    os.remove(outfile)
    return seqs

# -----------------------------------------------------------------

def GetHomologs(seq) :
    os.environ["BLASTDB"] = BLASTDB
    prog = config.BLASTDIR + "blastall"
    tmp = tempfile.mktemp()
    infile = tmp + ".blastin"
    outfile = tmp + ".blastout"

    fp = open(infile, "w")
    fp.write(seq)
    fp.close()
    
    cmd = "%s -p blastp -d nr -i %s -o %s" % (prog, infile, outfile)
    os.system(cmd)

    # Now we have to parse this file.  This code was adapted from
    # http://www.bioinformatics.org/bradstuff/bp/tut/Tutorial003.html#toc12
    b_fp = open(outfile, "r")
    b_parser = NCBIStandalone.BlastParser()
    b_record = b_parser.parse(b_fp)
    
    E_VALUE_THRESH = 0.04

    seqs = []
    for alignment in b_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                ss = hsp.sbjct.replace("-", "")
                seqs.append((alignment.title, ss))
                print "SUBJ" + hsp.sbjct
                print "MATCH" + hsp.match
    os.system("rm %s*" % tmp)
    return seqs

# -----------------------------------------------------------------

# Also see http://hg.wustl.edu/info/blast-program.html
# To run wu-blast you also need to run ~bpg/wu-blast/setdb nr

def MakeDatabase(name, seqs) :
    if 0 :
        ginums = {}
        for s in seqs.keys() :
            if s.startswith("gi|") :
                s1 = s
                s = s.split('|')[1]
                if ginums.has_key(s) : raise s
                ginums[s] = s1

    fasta.WriteSequences(seqs, name)
    RunCmd("%sformatdb -i %s -p T -o T" % (config.BLASTDIR, name))
    
# -----------------------------------------------------------------

def MakeDatabaseFromFile(name) :
    RunCmd("%sformatdb -i %s -p T -o T" % (config.BLASTDIR, name))
    RunCmd("%ssetdb %s" % (config.WUBLASTDIR, name))
    
# -----------------------------------------------------------------

if __name__ == "__main__" :
    if sys.argv[1] == "psiblast" :
        infile = sys.argv[2]
        prefix = sys.argv[3]
        params = " ".join(sys.argv[4:])
        seq = fasta.ReadOneSequence(infile)
        fasta.WriteSequence(prefix + ".tmp", "__query__", seq)
        cmd = "/home/bpg/blast/blastpgp -i %s.tmp -j 5 -m 6 -d /home/bpg/blast_db/nr_5_30_2003 -b 10000 -o %s.out %s" % (prefix, prefix, params)
        print cmd
        os.system(cmd)
        cmd = "/home/dek/python-2.2/bin/python2.2 /home/jechan/work/scripts/valid/blastParser.py %s.out %s.a2m" % (prefix, prefix)
        print cmd
        os.system(cmd)
        #aln = LoadAlignment("__query__", prefix + ".out")
        #print "loaded %d sequences in alignment" % len(aln)
        #fasta.WriteSequences(aln, prefix + ".a2m")
