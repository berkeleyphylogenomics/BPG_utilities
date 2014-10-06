#!/usr/bin/python

import fnmatch
import hmac
import os
import re
import string
import sys

from bpg.common import blosum62
from bpg.common.legacy_utils import *

ResidueNumber = { '-' : -1, '.': -1 }
NumberResidue = {}
i = 0
for c in "ACDEFGHIKLMNPQRSTVWY" :
    ResidueNumber[c] = i
    ResidueNumber[c.lower()] = i
    NumberResidue[i] = c
    i = i + 1

ThreeLetterCodes = {
    "ALA" : "A",
    "ARG" : "R",
    "ASN" : "N",
    "ASP" : "D",
    "CYS" : "C",
    "GLA" : "U",
    "GLN" : "Q",
    "GLU" : "E",
    "GLY" : "G",
    "HIS" : "H",
    "ILE" : "I",
    "LEU" : "L",
    "LYS" : "K",
    "MET" : "M",
    "ORN" : "O",
    "PHE" : "F",
    "PRO" : "P",
    "SER" : "S",
    "THR" : "T",
    "TRP" : "W",
    "TYR" : "Y",
    "VAL" : "V",
}

CodonUsage = {
    "TTT" : "F",
    "TCT" : "S",
    "TAT" : "Y",
    "TGT" : "C",
    "TTC" : "F",
    "TCC" : "S",
    "TAC" : "Y",
    "TGC" : "C",
    "TTA" : "L",
    "TCA" : "S",
    "TAA" : "*",
    "TGA" : "*",
    "TTG" : "L",
    "TCG" : "S",
    "TAG" : "*",
    "TGG" : "W",
    "CTT" : "L",
    "CCT" : "P",
    "CAT" : "H",
    "CGT" : "R",
    "CTC" : "L",
    "CCC" : "P",
    "CAC" : "H",
    "CGC" : "R",
    "CTA" : "L",
    "CCA" : "P",
    "CAA" : "Q",
    "CGA" : "R",
    "CTG" : "L",
    "CCG" : "P",
    "CAG" : "Q",
    "CGG" : "R",
    "ATT" : "I",
    "ACT" : "T",
    "AAT" : "N",
    "AGT" : "S",
    "ATC" : "I",
    "ACC" : "T",
    "AAC" : "N",
    "AGC" : "S",
    "ATA" : "I",
    "ACA" : "T",
    "AAA" : "K",
    "AGA" : "R",
    "ATG" : "M",
    "ACG" : "T",
    "AAG" : "K",
    "AGG" : "R",
    "GTT" : "V",
    "GCT" : "A",
    "GAT" : "D",
    "GGT" : "G",
    "GTC" : "V",
    "GCC" : "A",
    "GAC" : "D",
    "GGC" : "G",
    "GTA" : "V",
    "GCA" : "A",
    "GAA" : "E",
    "GGA" : "G",
    "GTG" : "V",
    "GCG" : "A",
    "GAG" : "E",
    "GGG" : "G",
}

# -----------------------------------------------------------------

def ReadSubfams(file) :
    fp = open(file, "r")
    seqs = {}
    subfams = {}
    subfam = ""
    while 1 :
        line = fp.readline()
        if not line :
            break
        
        if line.startswith("%subfamily ") :
            subfam = line[11:-1]
            subfams[subfam] = []
        elif line[0] == ">" :
            name = line[1:-1]
            seqs[name] = ""
            if subfam :
                subfams[subfam].append(name)
        else :
            seqs[name] = seqs[name] + line[0:-1]

    fp.close()
    return seqs, subfams

# -----------------------------------------------------------------

def ReadSequences(file, maxnum = 0, pats = [], notpats = []) :
    seqs = {}
    for pair in ReadSequencesList(file, maxnum, pats, notpats) :
        seqs[pair[0]] = pair[1]
    return seqs

# -----------------------------------------------------------------

def ReadOneSequence(file) :
    seqs = ReadSequencesList(file, 1)
    if not seqs : return None
    return seqs[0][1]

# -----------------------------------------------------------------

# Break up the multiple parts in a given name and return multiple components

def SimplifyNames(seqs) :
    nseqs = {}
    for s in seqs :
        ents = s.split('|')
        i = 0
        while i < len(ents) - 1 :
            if ents[i] in ("gi", "pdb", "sp", "ref", "lcl", "emb") :
                nseqs[ents[i] + '|' + ents[i + 1]] = s
                i = i + 1
            i = i + 1
    return nseqs

# -----------------------------------------------------------------

def SimplifiedNameOverlap(seqs1, seqs2, name1) :
    for s1 in seqs1 :
        if seqs1[s1] == name1 :
            if seqs2.has_key(s1) : return 1
    return 0

# -----------------------------------------------------------------

def OkToAppend(name, pats, notpats) :
    ok = 1
    if pats :
        ok = 0
        for p in pats :
            if fnmatch.fnmatch(name, p) :
                ok = 1
                break
    for p in notpats :
        if fnmatch.fnmatch(name, p) :
                ok = 0
                break
        
    return ok

# -----------------------------------------------------------------

def ReadSequencesList(file, maxnum = 0, pats = [], notpats = [], nop1 = 0) :
    fp = open(file, "r")
    seqs = []
    cur = ""
    name = "target"
    count = 0
    while 1 :
        line = fp.readline()
        if not line :
            if cur and OkToAppend(name, pats, notpats) :
                seqs.append((name, cur))
            break
        if line[0] == ">" :
            if cur and OkToAppend(name, pats, notpats) :
                seqs.append((name, cur))
                count = count + 1
                if maxnum and count >= maxnum : break
            name = line[1:-1].strip().replace('\t', ' ')
            if nop1 and name.startswith("P1;") : name = name[3:]
            cur = ""
        else :
            sub = line.rstrip().replace(' ', '').replace('\t', '').replace('*', '').replace('\r', '')
            cur = cur + sub
    fp.close()
    return seqs

# -----------------------------------------------------------------

def ReadFullNames(sfile) :
    names = {}
    cur = ""
    name = "target"
    count = 0
    for line in file(sfile) :
        if line[0] == ">" :
            line = line[1:-1]
            name = line.split()[0]
            names[name] = line
    return names

# -----------------------------------------------------------------

def ParseSeqs(string) :
    lines = string.split("\n")
    seqs = {}
    for line in lines :
        line = line.strip()
        if not line : continue
        if line[0] == ">" :
            name = line[1:]
            seqs[name] = ""
        else :
            line = line.replace(" ", "")
            line = line.replace("\t", "")
            seqs[name] = seqs[name] + line
    return seqs

# -----------------------------------------------------------------

def Consensus(seqs, subfams, lc_if_less = 0.0) :
    cons = {}
    for sf in subfams.keys() :
        ll = len(seqs[subfams[sf][0]])
        for ss in subfams[sf] :
            if ll != len(seqs[ss]) :
                raise "length mismatch: %s = %d, want %d" % (ss, len(seqs[ss]), ll)
            
        cs = ""
        for pos in range(0, ll) :
            count = {}
            totc = 0
            for n in subfams[sf] :
                c = seqs[n][pos]
                if c.islower() or c == '-' or c == '.' : continue
                    
                if count.has_key(c) :
                    count[c] = count[c] + 1
                else :
                    count[c] = 1
                totc = totc + 1
            bestc = -1
            best = '-'
            for c in count.keys() :
                if count[c] > bestc :
                    bestc = count[c]
                    best = c
            if totc and float(bestc) / totc < lc_if_less :
                best = best.lower()
            cs = cs + best
        cons[sf] = cs
    return cons

# -----------------------------------------------------------------

def ConsensusUnaligned(orig_seqs, subfams) :
    cons = {}
    seqs = {}
    for ss in orig_seqs.keys() :
        new = ""
        for c in orig_seqs[ss] :
            if c.islower() or c == '.' : continue
            new = new + c
        seqs[ss] = new
        
    for sf in subfams.keys() :
        ll = len(seqs[subfams[sf][0]])
        for ss in subfams[sf] :
            if ll != len(seqs[ss]) :
                raise "length mismatch: %s = %d, want %d" % (ss, len(seqs[ss]), ll)
            
        cs = ""
        for pos in range(0, ll) :
            count = {}
            matched = 0
            for n in subfams[sf] :
                c = seqs[n][pos]
                if c == '-' : continue
                
                if count.has_key(c) :
                    count[c] = count[c] + 1
                else :
                    count[c] = 1
                matched = matched + 1
                
            bestc = -1
            best = '-'
            for c in count.keys() :
                if count[c] > bestc :
                    bestc = count[c]
                    best = c
            if bestc > matched * 0.70 :
                cs = cs + best
            else :
                cs = cs + best.lower()
                
        cons[sf] = cs
        
    return cons

# -----------------------------------------------------------------
# Calculate consensus -- show most-frequent character in each position.
# If most-frequent character is in 85% of sequences, show in uppercase.
# Ignore lowercase and ".", but include "-" as a "21st character".

trivial_translation = string.maketrans('', '')
dotlowercase = '.' + string.lowercase
def ConsensusOfSequences( seqs, consensus_criterion=0.70, \
                                                        gap_as_21st_char_f=0 ) :
    if not seqs : 
        return []
    
    # Set of sequences without inserts.

    newseqs = []
    for n in seqs :
        s = seqs[n]
        new = s.translate(trivial_translation, dotlowercase)
        newseqs.append(new)

    ll = len(newseqs[0])
    for ss in newseqs :
        if ll != len(ss) :
            raise "length mismatch: %s = %d, want %d" % (ss, len(ss), ll)
        
    cs = ""
    matched = len( newseqs )
    cutoff = float( matched )*float( consensus_criterion ) 
    for pos in range(0, ll) :
        count = {}
        ###matched = 0
        for ss in newseqs :
            c = ss[pos]
                
            # When skipping gaps, need to know number matched.
            # REVISED: criterion for uppercase is all sequences, not just those
            # without gaps in this position.

            if not gap_as_21st_char_f and c == '-' : 
                continue

            ###matched += 1

            if count.has_key(c) :
                count[c] += 1
            else :
                count[c] = 1

            
        # If column is all gaps, will get "-" (whether or not skipping gaps).

        best = '-'
        bestc = -1
        for c in count.keys() :
            if count[c] > bestc :
                bestc = count[c]
                best = c

        if float( bestc ) > cutoff:
            cs = cs + best
        else :
            cs = cs + best.lower()

    return cs
  

# -----------------------------------------------------------------

def OnlyMatchesMultiple(seqs) :
    nseqs = {}
    ll = 0
    for key in seqs.keys() :
        mm = ""
        for c in seqs[key] :
            if c.isupper() or c == '-' :
                mm = mm + c
        nseqs[key] = mm
        if ll :
            if len(mm) != ll :
                raise "lengths don't match: %d vs %d" % (len(mm), ll)
        else :
            ll = len(mm)
    return nseqs

# -----------------------------------------------------------------

def OnlyResidues(seqs) :
    nseqs = {}
    for key in seqs.keys() :
        nseqs[key] = OnlySeqResidues(seqs[key])
    return nseqs

# -----------------------------------------------------------------

def OnlySeqResidues(str) :
    mm = ""
    for c in str :
        if c.isupper():
            mm = mm + c
        elif c.islower():
            mm = mm + c.upper()
    return mm

# -----------------------------------------------------------------
matchupperdash = re.compile("[A-Z\-]")   # compile for all Uppercase and Dashes (match states & dashes)
matchdash = re.compile("[\-]")
matchupper = re.compile("[A-Z]")
matchseq = re.compile("[A-Z|a-z]")

def OnlyMatchesOrDashes(seq) :
    return "".join(matchupperdash.findall(seq))  # return only match states and dashes

# -----------------------------------------------------------------

def LongestInsert(seq) :
    il = 0
    long = 0
    for s in seq :
        if s.islower() :
            il = il + 1
        else :
            if il > long : long = il
            il = 0
    if il > long : long = il
    return long

# -----------------------------------------------------------------

def OnlyMatchesOrDashesWithPositions(seq) :
    nseq = ""
    positions = []
    i = 1
    for c in seq :
        if c.isupper() :
            nseq = nseq + c
            positions.append(i)
            i = i + 1
        elif c == '-':
            nseq = nseq + c
            positions.append(0)
        elif c.islower() :
            i = i + 1
    return nseq, positions

# -----------------------------------------------------------------

def AlignmentLength(aseqs) :
    ll = 0
    for s in aseqs :
        l2 = len(OnlyMatchesOrDashes(aseqs[s]))
        if not ll :
            ll = l2
        elif ll != l2 :
            Error("length mismatch in alignment")
    return ll

# -----------------------------------------------------------------

def OnlyMatches(seq) :
##     nseq = ""
##     for c in seq :
##         if c.isupper():
##             nseq = nseq + c
##     return nseq
    return "".join(matchupper.findall(seq))  # return only match states and dashes

# -----------------------------------------------------------------

def FractionMatched(seqs) :
    nres = 0
    nmat = 0
    nmatdel = 0
    for key in seqs.keys() :
        for c in seqs[key] :
            if c.isupper():
                nmat = nmat + 1
                nmatdel = nmatdel + 1
            elif c == '-' :
                nmatdel = nmatdel + 1
            nres = nres + 1
    return nmat / float(nres), nmatdel / float(nres), 

# -----------------------------------------------------------------

def WriteSequences(seqs, file=None) :
    str = ""
    if file:
        fp = open(file, "w")
        if type(seqs) == type([]) :
            order = []
            for name, val in seqs :
                fp.write(">%s\n%s\n" % (name, val))
        else :
            order = seqs.keys()
            for key in order :
                fp.write(">%s\n%s\n" % (key, seqs[key]))

        fp.close()
        return order

    else:
        if type(seqs) == type([]) :
            order = []
            for name, val in seqs :
                #str += ">%s\n%s\n" % (name, val)  
                str = "".join("%s>%s\n%s\n" % (str, name, val)) # faster if str >= 10K
        else :
            order = seqs.keys()
            for key in seqs.keys() :
                #str += ">%s\n%s\n" % (key, seqs[key])
                str = "".join("%s>%s\n%s\n" % (str, key, seqs[key]))
        return str

# -----------------------------------------------------------------

def WriteSequencesPir(seqs, annotations, file) :
    fp = open(file, "w")
    for key in seqs.keys() :
        fp.write(">P1;%s\n%s\n%s*\n" % (key, annotations[key], seqs[key]))
    fp.close()

# -----------------------------------------------------------------

def WriteSequence(file, name, val) :
    fp = open(file, "w")
    fp.write(">%s\n" % (name))
    for i in range(0, len(val), 72) :
        fp.write("%s\n" % val[i:i+72])
    fp.close()

# -----------------------------------------------------------------

def RemapNames(seqs) :
    i = 0
    map = {}
    nseqs = {}
    for s in seqs.keys() :
        s2 = "tmp|%d" % i
        nseqs[s2] = seqs[s]
        map[s2] = s
        map[s] = s2
        i = i + 1
    return nseqs, map

# -----------------------------------------------------------------

def SeqFromPdbAtoms(ff, chain) :
    seq = ""
    pos = []
    fp = file(ff)
    for line in fp :
        if not line.startswith("ATOM") : continue
        res = line.split()[3]
        if chain :
            ch = line.split()[4]
            if ch != chain : continue
            p = int(line.split()[5])
        else :
            p = int(line.split()[4])
        if not pos or p != pos[-1] :
            pos.append(p)
            seq = seq + ThreeLetterCodes[res]
    return seq, pos
    
# -----------------------------------------------------------------

def SeqFromPdbSeqres(ff, chain) :
    seq = ""
    pos = []
    stuff = ""
    fp = file(ff)
    for line in fp :
        if not line.startswith("SEQRES") : continue
        if chain and chain != line.split()[2] : continue
        stuff = stuff + line
    return ThreeLetterToSeq(stuff)

# -----------------------------------------------------------------

def ThreeLetterToSeq(text) :
    seq = ""
    stuff = text.split()
    for word in stuff :
        word = word.upper()
        if ThreeLetterCodes.has_key(word) :
            seq = seq + ThreeLetterCodes[word]
    return seq

# -----------------------------------------------------------------

def TrimEnds(seq) :
    nseq = ""
    insert = ""
    for c in seq :
        if c.isupper() :
            if nseq : nseq = nseq + insert
            insert = ""
            nseq = nseq + c
        elif c.islower() :
            insert = insert + c.upper()
    return nseq
      
# -----------------------------------------------------------------

# These probably shouldn't be different functions

def TrimEnds2(seq) :
    nseq = ""
    insert = ""
    for c in seq :
        if c.isupper() or c == '-' :
            if nseq : nseq = nseq + insert
            insert = ""
            nseq = nseq + c
        elif c.islower() :
            insert = insert + c
    return nseq
      
# -----------------------------------------------------------------
#This function will trim periods ('.'), dashes ('-'), and lowercase characters
#from the end and from the beginning until an uppercase character is hit
def TrimEnds3(seq) :
    nseq = ""
    insert = ""
    for c in seq :
        if c.isupper():
            if nseq :
                nseq = nseq + insert + c
                insert = ""
            else:
                insert = ""
                nseq = c
        else:
            insert = insert + c
    return nseq
      
# -----------------------------------------------------------------

def PairwiseIdentity(seq1, seq2) :
    seq1 = OnlyMatchesOrDashes(seq1)
    seq2 = OnlyMatchesOrDashes(seq2)
    if len(seq1) != len(seq2) : Error("string lengths don't match")
    id = 0
    for i in range(len(seq1)) :
        if seq1[i] == seq2[i] and seq1[i] != '-' : id = id + 1
    return float(id) / len(seq1)

# -----------------------------------------------------------------

def DNAtoProtein(string) :
    dna = ""
    for word in string.split() :
        if not word.isalpha() : continue
        dna = dna + word
    seq = ""
    
    for i in range(0, len(dna), 3) :
        codon = dna[i:i+3].upper().replace('U', 'T')
        if CodonUsage.has_key(codon) :
            r = CodonUsage[codon]
        else :
            r = "?"
        seq = seq + r
    return seq

# -----------------------------------------------------------------

# Create six amino acid sequences from one DNA sequence -- based on 3 forward 
# and 3 reverse translations.

def DNAtoProtein6( seqfile ) :
    EST_sequences = ReadSequencesList( seqfile )
    sequences = []
    for EST_seq in EST_sequences :
        dna = EST_seq[1]
        rdna = ""
        for c in dna : 
            rdna = c + rdna

        for i in [0, 1, 2] :
	    header = "FORWARD%d " % i + EST_seq[0]
            res = DNAtoProtein(dna[i:]).replace("?", "")
	    sequences.append( ( header, res ) )

   	    header = "REVERSE%d " % i + EST_seq[0]
            rres = DNAtoProtein(rdna[i:]).replace("?", "")
   	    sequences.append( ( header, rres ) )

    return sequences

# -----------------------------------------------------------------

def CheckDnaMatchesProtein(proseq, dnaseq) :
    for i in range(len(proseq)) :
        res = CodonUsage.get(dnaseq[i:i+3].upper(), "?")
        if res != proseq[i].upper() : return 0
    return 1
    
# ----------------------------------------------------------------

def GetMatches(aseq, cseq) :
    if len(aseq) != len(cseq) :
        print aseq
        print cseq
        print "length doesn't match"
    match = ""
    for i in range(0, min(len(aseq), len(cseq))) :
        a = aseq[i].upper()
        c = cseq[i].upper()
        if a == c and a != ' ' and a != '-' :
            m = c
        else :
            v = blosum62.Value(a, c)
            if v > 0 :
                m = '+'
            else :
                m = ' '
        match = match + m
    return match

# -----------------------------------------------------------------

def CopyAndAddTag(fromfile, tofile, longnames, pref, vals) :
    seqs = ReadSequences(fromfile)
    nseqs = {}
    for s in seqs :
        # This might be slow if longnames has lots of entries ....
        s2 = s
        for ln in longnames :
            if ln.find(s2) >= 0 :
                s2 = ln
                break
        if vals.has_key(s) :
            #spl = s2.split('|')
            #spl = spl[0:-1] + [pref, "%s" % vals[s]] + spl[-1:]
            #nseqs[string.join(spl, '|')] = seqs[s]
            s2 = s2 + "|evalue|%s" % vals[s]
        nseqs[s2] = seqs[s]
    WriteSequences(nseqs, tofile)
    print "writing", tofile
    
# -----------------------------------------------------------------

def GetTag(name, tag) :
    next = 0
    for item in name.split('|') :
        if next : return item
        if item == tag : next = 1
        spl = item.split('=')
        if len(spl) > 1 and spl[0] == tag :
            return spl[1]
    return None

# -----------------------------------------------------------------

def NumberLines(seq) :
    n = 1
    nc = 0
    for s in seq :
        if s != ' ' and s != '.'  and s != '-' : nc = nc + 1
    num1 = ""
    num2 = ""
    nd = 0
    for s in seq :
        if s == ' ' or s == '.'  or s == '-' :
            if not nd :
                num1 = num1 + ' '
            else :
                nd = nd - 1
            num2 = num2 + ' '
        else :
            if n == 1 or n == nc or (n % 10 == 0 and n < nc - 5) :
                nstr = "%d" % n
                nd = len(nstr) - 1
                num1 = num1 + nstr
                num2 = num2 + '|'
            else :
                if not nd :
                    num1 = num1 + ' '
                else :
                    nd = nd - 1
                num2 = num2 + ' '
                    
            n = n + 1
    return num1, num2

# -----------------------------------------------------------------

def ExpandAlignment(raw_cseq, aseq) :

    if len(OnlyMatchesOrDashes(aseq)) != len(raw_cseq) :
        raise "length mismatch"
    
    # Now expand the consensus sequence so that it has gaps for the lower
    # case positions in the aligned sequence.
    cseq = ""
    orig_raw_cseq = raw_cseq
    for a in aseq :
        if a.islower() :
            cseq = cseq + " "
        else :
            if raw_cseq == "" :
                print "ASEQ is %s" % aseq
                raise "ERROR: raw_cseq == blank, orig = %s" % (orig_raw_cseq)
                cseq = cseq + " "
            else :
                cseq = cseq + raw_cseq[0]
                raw_cseq = raw_cseq[1:]

    if raw_cseq :
        raise "Error: left over stuff: %s" % raw_cseq
        #sys.exit(1)
        
    if len(aseq) != len(cseq) :
        print "Error: mismatched lengths:\n%s\nvs\n%s" % (aseq, target)
        #sys.exit(1)
    return cseq

# -----------------------------------------------------------------

def PrettyAlign(seqlist, consensuslist):

    returnedseqs = {}
    for seq in seqlist:
        seqlist[seq] = list(seqlist[seq].replace('.',''))
        returnedseqs[seq] = ""
        
    for seq in consensuslist:
        consensuslist[seq] = list(consensuslist[seq].replace('.',''))
        returnedseqs[seq] = ""


    notdone = 1
    while notdone:

        insert = 0
        for seq in seqlist:
            if seqlist[seq] != [] and seqlist[seq][0].islower():
                insert=1
                break;

        
        notdone = 0
        for (seqs, consensus) in [(consensuslist, 1), (seqlist, 0)]:

            for seq in seqs:

                remaining = seqs[seq]
                output = returnedseqs[seq]
                
                if remaining == []:
                    returnedseqs[seq] = output + '.'
                    continue
                else:
                    notdone = 1
                
                if not insert or (remaining[0].islower() and not consensus):
                    seqs[seq] = remaining[1:]
                    returnedseqs[seq] = output + remaining[0]
                else:
                    returnedseqs[seq] = output + '.'

    for seq in returnedseqs:
        returnedseqs[seq] = returnedseqs[seq][:-1]
    return returnedseqs


# ------------------------------------------------------------------------------
def is_dna_or_rna(seq):
    """Verifies wether a seq is DNA or RNA (composed only of [acgtuxn] chars).
       Inputs:
       seq = query sequence.
       Outputs:
       is_dna_rna_f = boolean indicating whether sequence is dna or rna.
    """
    is_dna_rna_f = False
    if re.compile('[^a-zA-Z]+').match(seq) :
        print seq
        raise ValueError, 'Sequence argument <seq> has invalid characters.'
    if len(seq) < 1 :
        raise ValueError, 'Invalid empty sequence!'
    if re.compile('^[acgtuxnACGTUXN]+$').match(seq):
        is_dna_rna_f = True
    else:
        is_dna_rna_f = False
    return is_dna_rna_f

# ------------------------------------------------------------------------------
def is_mostly_dna(seq, ratio):
    """Verifies wether a seq is mostly DNA or RNA (composed mostly of [acgtu] chars).
       Inputs:
       seq = query sequence.
       ratio = min ratio between dna_content and seq_length that will cause seq
               to be flagged as dna.
       Outputs:
       is_mostly_dna_f = boolean indicating whether sequence is mostly dna or rna.
    """
    is_mostly_dna_f = False
    if ratio < 0 or ratio > 1 :
        print 'ratio = %d' % ratio
        raise ValueError, '<ratio> has to be between 0 and 1.'
    if re.compile('[^a-zA-Z]+').match(seq) :
        print seq
        raise ValueError, 'Sequence argument <seq> has invalid characters.'
    if len(seq) < 1 :
        raise ValueError, 'Invalid empty sequence!'
    dna_letters = 0
    for letter in seq :
        if re.compile('[acgtuACGTU]').match(letter) :
            dna_letters = dna_letters + 1
    final_ratio = float(dna_letters)/float(len(seq))
    if final_ratio > ratio :
        is_mostly_dna_f = True
    return is_mostly_dna_f

# -----------------------------------------------------------------

if __name__ == "__main__" :
    if sys.argv[1] == "fixnames" :
        seqs = ReadSequencesList(sys.argv[2])
        seqs2 = []
        for s in seqs :
            nname = s[0]
            n = nname.split('|')
            if len(n) == 1 :
                n = n[0].split(' ')
            if len(n) > 1 :
                nlist = n[0:2] + [hmac.new(s[0]).hexdigest()]
                nname = "_".join(nlist)
            seqs2.append((nname, s[1]))
        WriteSequences(seqs2, sys.argv[3])

    elif sys.argv[1] == "copy" :
        seqs = ReadSequencesList(sys.argv[2])
        WriteSequences(seqs, sys.argv[3])
