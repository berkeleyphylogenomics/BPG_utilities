#! /usr/bin/env python
# alignments_lib.py       Terry Farrah    March 2008
#
# Library of classes and functions useful for the
# handling alignments within the SHMM-SHMM alignment accuracy project.

import pdb
import sys, time
import math, random
import heapq
import os, glob, re, commands, string
import BPG_common.fasta
from matchmaker.shmm_shmm_lib import *

_GAP = sys.maxint
_gapchar = '-'
_dotchar = '.'
_gapchars = set([_gapchar, _dotchar])

class AlignmentLibError(Exception):
   pass

#=======================================================================
# Alignment data structures and functions

class Alignment:
   """A pairwise algnmnt, represented by 2 strings of equal length with gaps.

   """

   def __init__(self):
      self.len = 0   # the length of the alignment (not nec. of each seq)
      self.fname = None   # filename of alignment
      self.dirname = None   # pathname of directory containing alignment
      self.pairsNoGaps = set()   # set of aligned residues (i,j);
                           # indexing starts at 0, (i,j) is a tuple
      self.pairsInclGaps = set() # .pairsNoGaps plus indel pairs
      self.cs_score = None  # Cline shift score relative to reference
      self.reference = None # reference alignment
      # the fields below are pairs, one for each sequence in the alignment
      self.name = ["",""]   # names of sequences
      self.line = ["",""]    # lines of alignment
      self.seq = ["",""]  # raw sequences (with gap characters removed)
      self.seqlen = ["",""]   # lengths of sequences
      # for each seq k={0,1}, seqpairs[k][i]=j if residue i is paired
      # with residue j in other seq, else =_GAP if unpaired.
      self.seqpairs = [[],[]] 

      # the fields below are used during clustering
      self.cluster = None   # cluster this alignment currently belongs to
      self.alignmentPairs = []  # list of alignmentPair data structures
                                # representing pairs involving this algnmnt

   def readFromFile(self, filename, pair=None):
      """read a pairwise alignment file into an Alignment data structure

      """
      """
      # 7/15/08: replaced this function with one that can read either compressed
      # or uncompressed files
      alignmentList = BPG_common.fasta.ReadSequencesList(filename)
      if len(alignmentList) != 2:
         raise AlignmentLibError, \
          "Did not read exactly 2 fasta sequences from file %s" \
           % filename
      """
      alignmentList = read_alignment_file(filename)
      alignmentLine = [ alignmentList[0][0], alignmentList[0][1],
         alignmentList[1][0], alignmentList[1][1] ]
      self.readFromStrings(alignmentLine, filename)
      # try to get the Cline shift score from separate file in same dir
      # strip off final "original.afa" and add "CS.out"
      score_filename = filename[:-12] + "CS.out"
      try:
         sfile = open(score_filename, "r")
         scoreline = sfile.readlines()[-1] #get last line of file
         tokens = scoreline.rsplit("CS=")
         cs_score_token = tokens[1].split(";")[0]
         self.cs_score = float(cs_score_token)
      except: pass


   def readFromStrings(self, alignmentLine, filename=None):
      ''' read fasta format pairwise alnmt from 4 strings and store in self.

      '''
      if filename != None:
         self.fname = os.path.basename(filename)
         self.dirname = os.path.dirname(filename)

      for i in [0,1]:
         self.name[i] = alignmentLine[i*2].strip()
         self.line[i] = alignmentLine[i*2+1].strip()

      for i in [0,1]:
         self.seq[i] = \
            self.line[i].replace(_gapchar,"").replace(_dotchar,"").upper()
         seqlen = len(self.seq[i])
         if seqlen >= _GAP:
            raise AlignmentLibError, \
                 "Sequence %d (%s) in %s longer than %d, definition of gap." \
                  %  (i+1, self.name[i], filename, _GAP)
         self.seqlen[i] = seqlen
         for k in range (self.seqlen[i]):
            self.seqpairs[i].append(None)

      if len(self.line[0]) != len(self.line[1]):
         raise AlignmentLibError, \
           "Alignment.readFromFile(): algnmt line lengths differ \
for seqs %s, %s  in file %s" \
             % (self.name[0], self.name[1], filename)
      self.len = len(self.line[1])
      self.positionPairsFromLines()


   def equivalent(self,a):
      ''' Test to see if two alignments are equivalent.

      Alignments are equivalent if (a) they are for the same seqs,
      in the same order, and (b) they have the same set of aligned pairs.
      '''
      return (self.hasSameSeqsAs(a) and self.pairsInclGaps == a.pairsInclGaps)
      

   def hasSameSeqsAs(self, a):
      '''Return True if a is an alignment for the same seqs, in same order.

      Assumes self.seq has been filled in.
      '''
      return (self.seq[0]==a.seq[0] and self.seq[1]==a.seq[1])


   def positionPairsFromLines(self):
      """ List my aligned pairs i,j (starting at 0) and store in self

      Assumes that self.seqpairs has been created, but is empty.
      Also create seqpairs list.
      """
      i = j = 0
      for k in range(self.len):
         if (self.line[0][k] not in _gapchars):
            if (self.line[1][k] not in _gapchars):
               self.pairsInclGaps.add((i,j))
               self.pairsNoGaps.add((i,j))
               self.seqpairs[0][i] = j
               self.seqpairs[1][j] = i
               j = j+1
               i = i+1
            else:
               self.pairsInclGaps.add((i,_GAP))
               self.seqpairs[0][i] = _GAP
               i = i+1
         elif (self.line[1][k] not in _gapchars):
            self.pairsInclGaps.add((_GAP,j))
            self.seqpairs[1][j] = _GAP
            j = j+1

   def seqPairsFromLines(self):
      """ From alignment lines, create self.seqpairs.

      Assumes that self.seqpairs has NOT been created.
      positionPairsFromLines does same, and also
      creates self.pairsInclGaps and self.pairsNoGaps.
      These two functions could be consolidated.
      """
      for i in [0,1]:
         for k in range (self.seqlen[i]):
            self.seqpairs[i].append(None)
      i = j = 0
      for k in range(self.len):
         if (self.line[0][k] not in _gapchars):
            if (self.line[1][k] not in _gapchars):
               self.seqpairs[0][i] = j
               self.seqpairs[1][j] = i
               j = j+1
               i = i+1
            else:
               self.seqpairs[0][i] = _GAP
               i = i+1
         elif (self.line[1][k] not in _gapchars):
            self.seqpairs[1][j] = _GAP
            j = j+1

   def alignmentLinesFromPairs(self, verbose=False):
      '''Create standard alignment display from set of pairs.

      Assumes self.pairsInclGaps and self.seq[] are filled in.
      "Standard display" is two strings with residue and gap characters.
      Store in self.line[].
      Inverse of positionPairsFromLines()
      12/20/07: currently breaks up unaligned regions -- not incorrect,
        but undesirable.
      '''

      # pick a pair from the set that includes the next position
      # in either sequence, until done.
      # (error messages below indicate program bug, not problem with data.)

      # I tried instead copying the set into a set, but it said
      # I couldn't iterate over a non-sequence.
      pairList = list(self.pairsInclGaps)

      # I think the code below could be made more concise

      i = j = 0               # index of next residue in each seq
      if verbose:
         print "Ordered aligned pairs:\n"
      while i < self.seqlen[0] or j < self.seqlen[1]:
         # brute force search for the pair we want
         foundNext = False
         for pair in pairList:
            if pair[0]==i:
               if pair[1]==j:
                  # found a pair with a residue in each position
                  self.line[0] = self.line[0] + self.seq[0][i]
                  self.line[1] = self.line[1] + self.seq[1][j]
                  i = i+1
                  j = j+1
                  pairList.remove(pair)
                  foundNext = True
                  if verbose:
                     print pair
                  break
               elif pair[1] == _GAP:
                  # found a pair with a gap in second position
                  self.line[0] = self.line[0] + self.seq[0][i]
                  self.line[1] = self.line[1] + _gapchar
                  i = i+1
                  pairList.remove(pair)
                  foundNext = True
                  if verbose:
                     print pair
                  break
            if pair[1]==j:
               if pair[0] == _GAP:
                  # found a pair with a gap in first position
                  self.line[0] = self.line[0] + _gapchar
                  self.line[1] = self.line[1] + self.seq[1][j]
                  j = j+1
                  pairList.remove(pair)
                  foundNext = True
                  if verbose:
                     print pair
                  break
         if not foundNext:
            print >> sys.stderr, \
       "alignmentLinesFromPairs: didn't find pair for i=%d or j=%d, file %s" \
                % (i,j, self.fname)
            self.line[0] = self.line[0] + '?'
            self.line[1] = self.line[1] + '?'
            i = i + 1
            j = j + 1
      if len(pairList) > 0:
         print >> sys.stderr, \
                  "pairList not empty at end of alignmentLinesFromPairs()."
         print >> sys.stderr, pairList
      self.len = len(self.line[0])
      return()

   def pairsNeededForCoverage(self, coverage):
      '''return # of pairs needed to cover given fraction of shorter seq'''
      #if self.seqlen[0] < self.seqlen[1]:
         #return int(coverage*self.seqlen[0])
      #else: return int(coverage*self.seqlen[1])
      return (10)


   def printLines(self, file=sys.stdout):
      '''print alignment in standard format '''
      print >> file, self.line[0]
      print >> file, self.line[1]

   def printPIRforModeller(self, file=sys.stdout, reverseOrder=False):
      '''print alignment in PIR format ready for input to Modeller

      3/15/08: written today but NOT TESTED; MAY NOT BE NEEDED.
      first seq will be template, second will be target.
      if reverseOrder is True, then reverse
      '''
      print >> file, "C;"
      print >> file, ""
      for i in [0,1]:
         print >> file, ">p1;%s" & self.name[i]
         if (reverseOrder and i==1) or (not reverseOrder and i==0):
            print >> file, "structure:%s" % self.name[i]
         else:
            print >> file, "sequence"
         print >> file, "%s*" % self.seq[i]

           

class AlignmentPlotPoint:
   """A point on an alignment plot.

   Each point contains the following information:
   data: 0=nothing here,  1=align, 2=upgap, 3=downgap
   pos: position in first sequence, in second sequence
   alignments:
     list of names of alignments that contribute to an "align" point.
   seqSegment: shows a portion of a sample alignment that spans this point
   """

   def __init__(self):
      self.data = 0
      self.pos = [None,None]
      self.alignments = []
      self.seqSegment = ["",""]

def checkAlignmentFileName(f):
   '''if name has certain substrings, return True.

   f = file basename
   return True if file is an alignment file we want to use
   '''

   # must end in .afa
   if not f.endswith(".afa"):
      return(False)
   # return True if contains YL or hmm
   #return ( f.find("YL") != -1 or f.find("hmm")!= -1)


#=======================================================================
# Alignment metrics
# Used by hierarchical agglomerative clustering.
# Input: two alignments  Output: distance (number)

def distanceBetweenAlignments(a1, a2, method="1-Qcombined"):
   ''' Return dist (normalized to 0..1) btwn two algnmnts using method.

   '''
   if method=="1-Qcombined":
      return QcombinedMetric(a1, a2)
   elif method=="AMAmetric":
      return(normalizedAMAmetric(a1, a2))
   elif method=="ClineShift":
      return(normalizedClineShiftDistance(a1, a2))
   else:
      raise AlignmentLibException, "Unknown metric %s" % method


def QcombinedMetric(a1, a2):
   ''' Convert Q-combined to a metric by subtracting from one. '''
   return 1-Qcombined(a1, a2)

def Qcombined(a1, a2, countIndels=False):
   ''' Return the Yona-Levitt Q-combined similarity for two alignments.

    This is the number of shared aligned pairs divided by the total
    number of aligned pairs in each alignment. Indels are NOT counted
    in Yona and Levitt's formulation.

    Yona, G. and Levitt, M. 2002. Within the twilight zone: A sensitive
    profile-profile comparison tool based on information theory.
    J. Mol. Biol. 315: 1257-1275.
   '''

   if countIndels:
      return(float(len(a1.pairsInclGaps.intersection(a2.pairsInclGaps))) /
            float(len(a1.pairsInclGaps.union(a2.pairsInclGaps))))
   else:
      lenPairUnion = len(a1.pairsNoGaps.union(a2.pairsNoGaps))
      if lenPairUnion == 0: return (0)  # avoid divide by zero
      return(float(len(a1.pairsNoGaps.intersection(a2.pairsNoGaps))) /
            float(lenPairUnion))


def normalizedAMAmetric(a1, a2):
   ''' Normalize AMA  metric to between 0 and 1.'''
   return(float(AMAmetric(a1,a2)) /
         float(len(a1.seq[0]) + len(a2.seq[0])))

def AMAmetric(a1, a2):
   ''' Return metric for two algnmts as in Schwartz, Myers, & Pachter

   This is the total number of residues from both seqs that
   are aligned differently in the two alignments.
   Schwartz, Meyers, and Pachter, "Alignment Metric Accuracy", 2006,
   formula (14):
          distance = len(s1) + len(s2) - 2(# shared pairs) -
                     # shared insertions - # shared deletions
   Could implement more cleanly using sets, as in Qcombined().
   '''

   nSharedPairs = nSharedIns = nSharedDel = 0
   # look at each position in the first sequence
   for i in range (a1.seqlen[0]):
      if a1.seqpairs[0][i] != _GAP:
         if a2.seqpairs[0][i] != _GAP:
            # if the 2 pairs are identical, we have a shared pair
            if a1.seqpairs[0][i] == a2.seqpairs[0][i]:
               nSharedPairs = nSharedPairs + 1
      # case where both alignments have a deletion here
      elif a2.seqpairs[0][i] == _GAP:
         nSharedDel = nSharedDel + 1
   # look at each position in the second sequence
   for i in range (a1.seqlen[1]):
      # both alignments have an insertion here
      if a1.seqpairs[1][i] == _GAP and a2.seqpairs[1][i] == _GAP:
         nSharedIns = nSharedIns + 1

   return (a1.seqlen[0] + a1.seqlen[1] - 2*nSharedPairs -
      nSharedIns - nSharedDel)


def normalizedClineShiftDistance(a1,a2,epsilon=0.2):
   '''Cline shift distance normalized to range 0..1'''
   return (clineShiftDistance(a1,a2,epsilon)/ (1+epsilon))


def clineShiftDistance(a1,a2,epsilon=0.2):
   '''1.0 minus Cline shift score. Not sure this is a metric.'''
   return ( 1.0 - clineShiftScore(a1,a2,epsilon))


def clineShiftScore(a1,a2,epsilon=0.2):
   '''Return shift similarity score as defined by M. Cline
   
   Shift score is the sum of score s(r)
   for all residues r that are paired in both alignments,
        s(r) = (1.2 / 1+shift(r)) -.2,
   divided by (# aligned pairs in a1 + # aligned pairs in a2)
   '''
   numerator = 0.0
   for k in [0,1]:                  # for each sequence
      for i in range(a1.seqlen[k]): # for each residue
         r1 = a1.seqpairs[k][i]
         r2 = a2.seqpairs[k][i]
         # if residue is paired in both alignments, incorporate into numerator
         if r1 != _GAP and r2 != _GAP:
             shift = abs(r1-r2)
             numerator = numerator + \
               (1 + epsilon) / float(1 + shift) - epsilon
   denominator = len(a1.pairsNoGaps) + len(a2.pairsNoGaps)
   score = numerator / denominator
   return(score)


def testMetrics(alignments):
   ''' Calculate/tally distances between each pair of alignments in a list.

   Purpose: discover the distribution of distances within the set,
   and compare the distributions for the different metrics.
   Jan. 2008. 
   '''

   # I know there must be a better way to initialize these arrays
   m = alignments[0].seqlen[0]
   n = alignments[0].seqlen[1]
   AMATally=[]
   QcombinedTally=[]
   for i in range(m+n):
      AMATally = AMATally + [0]
      QcombinedTally = QcombinedTally + [0]
   
   for i in range (nAlignments-1):
      for j in range (i+1, nAlignments):
   
         # make sure the sequences in these 2 alignments are identical.
         canCompare = True
         for k in [0,1]:
            if (alignments[i].seq[k] != alignments[j].seq[k]):
               print >> sys.stderr, "Sequences %d in alignments %d and %d differ" % (k+1, i, j)
               print >> sys.stderr, alignments[i].fname
               print >> sys.stderr, alignments[i].name[k]
               print >> sys.stderr, alignments[i].seq[k]
               print >> sys.stderr, alignments[j].fname
               print >> sys.stderr, alignments[j].name[k]
               print >> sys.stderr, alignments[j].seq[k]
               canCompare = False
   
         if canCompare:
            qmetric = 1 - Qcombined(alignments[i], alignments[j])
            # scale 1 - Qcombined so that it's an integer ranging from 0 to m+n, like AMA  metric
            qmetricInt = int(round(qmetric * float(m+n)))
            pmetric = AMAmetric(alignments[i], alignments[j])
            #  scale AMA metric so it ranges from 0 to 1, like 1 - Qcombined
            pmetricFloat = float(pmetric) / float(m+n)
   
            # tally the metrics to see their distribution (development purposes)
            AMATally[pmetric] = AMATally[pmetric] + 1
            QcombinedTally[qmetricInt] = QcombinedTally[qmetricInt] + 1

   # print tallies into a tab-separated table for import into Excel
   # (implemented approx. 12/21/07 to create a scatter plot to show
   #  distribution of distances & compare distributions for two
   #  metrics.)
   print "%s vs. %s" % (alignments[0].name[0], alignments[0].name[1])
   print "tally\tqmetricInt\tpmetric"
   for i in range(m+n):
      print "%d\t%d\t%d" % (i, AMATally[i], QcombinedTally[i])
   print nAlignments
