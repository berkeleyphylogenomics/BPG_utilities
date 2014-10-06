#! /usr/bin/env python
# coding: iso-8859-15

""" consensus_alignment.py       Terry Farrah    December 13, 2007

Input:
 - set of pairwise alignment files. Each alignment
   is between the same pair of sequences. 
Output:
 - list of consensus sequences representing alignment clusters

Method:

 Cluster alignments using:
  1)  heirarchical agglomerative clustering
      with (a) joining to centroid (slow; disabled 2/18/08) or
         (b) single linkage.
  2) heuristic hash on residue pairs
  3) (2) recursively down to moderate cluster size, then (1)

 Metric for hierarchical agglomeration is one of following:
  1) Q-combined (shared aligned pairs divided by total aligned pairs)
      This is a Jaccard dissimilarity.
  2) Lior Pachter's AMA alignment metric:
      total number of residues from both seqs that are aligned
      differently in the two alignments
  3) Jaccard, including gaps in the def. of "aligned pair"
        (not yet implemented)

 Form consensus for each cluster using one of following:
  1) greedy algorithm over set of all aligned pairs (including gaps)
       (fails on degenerate alignment sets)
  2) consensus = A* = minarg(median(A*,A')) where A' belongs to cluster
  3) AMAP annealing (not yet fully implemented)
"""

import sys, os, glob, pdb
import time
import Numeric
from matchmaker.alignments_lib import *
from matchmaker.cluster_consensus_lib import *
from matchmaker.shmm_shmm_lib import *
from optparse import OptionParser

_coverageRequirement = 0.35

class ClusteringAndConsensusError(AlignmentLibError):
   pass

#====================
# Process command line
#====================
parser = OptionParser()
parser.add_option("--seedpair", dest="seedpair",
      help="matchmaker pair of seeds of form seedx_seedY. " + \
           "Either --seedpair or --files mandatory.",
      metavar="PAIR")
parser.add_option("-f", "--files", dest="files",
      help="files containing input alignments; enclose in quotes if using *, ?, etc.",
      metavar="FILES")
parser.add_option("--nr", dest="use_nr", default=False, action="store_true",
      help="cluster nonredundant alignments" + \
           "stored in <pair>_align/<pair>_nr.csv " + \
            "(must use in combination with --seedpair)")
parser.add_option("-v", "--verbose", dest="verbose", default=False,
                  action="store_true",
      help="print extra status messages to stdout")
parser.add_option("-q", "--quiet", dest="quiet", default=False,
                  action="store_true",
      help="don't print anything except the result")
parser.add_option("-r", "--reference", dest="ref_alignment_file",
      metavar="REF",
      help="compare all consensus alignments to alignment in this file")
parser.add_option("-m", "--metric", dest="metric", default="1-Qcombined",
      metavar="METRIC",
      type="choice",
      choices=["1-Qcombined","AMAmetric"],
      help="Specify metric as 1-Qcombined (default) or AMAmetric")
parser.add_option("-c", "--clustering", dest="clusterMethod",
      metavar="METHOD",
      type="choice",
      choices=["singleLinkage","centroidJoining", "heuristicHash", "dualPhase"],
      default="dualPhase",
help="Specify clustering as singleLinkage, centroidJoining (slower), heuristicHash, or dualPhase (default; heuristic hash down to manageable-size clusters, followed by agglomerative within clusters)")
parser.add_option("-d", "--distanceThreshold", dest="distanceThreshold",
      metavar="DISTANCE",
      type="float",
      default=0.05,
      help = "Stop clustering when distances between clusters is greater than this; should be between 0 and 1; default 0.05")
parser.add_option("-x", "--maxClusters", dest="maxClusters",
      metavar="COUNT",
      default=500,
      type="int",
      help="Max # of clusters; supercedes distanceThreshold. Default=500.")
parser.add_option("-p", "--pairFrequency", dest="pairFrequency",
      metavar="FRACTION",
      type="float",
      default=0.6,
help="For heuristic hash: set approx. residue pair frequency to hash on (default 0.6).")
parser.add_option("-s", "--skip", dest="nPairsToSkip",
      metavar="N",
      type="int",
      default=4,
help="For heuristic hash: hash on every Nth pair in pair list (default 4)")
parser.add_option("-k", "--consensusMethod", dest="consensusMethod",
      metavar="METHOD",
      type="choice",
      choices=["greedy","minMedian","none"],
      default="minMedian",
      help = "Consensus method; choose minMedian (default), greedy (fails on diverse input), or none (in this case minMedian will be used during clustering)")

(options, args) = parser.parse_args()

# check that necessary options are given, and that values are valid
# assign option values to variables
if not (options.files or options.seedpair):
   parser.error("Either --seedpair or --files required")
if options.files and options.seedpair:
   parser.error("Only one of --seedpair or --files allowed")
if options.files and options.use_nr:
   parser.error("--nr can only be used in combination with --seedpair")
files = options.files
seedpair = options.seedpair
use_nr = options.use_nr
metric = options.metric
clusterMethod = options.clusterMethod
consensusMethod = options.consensusMethod
verbose = options.verbose
quiet = options.quiet
ref_alignment_file = options.ref_alignment_file
maxClusters = options.maxClusters
threshold = options.distanceThreshold
nPairsToSkip = options.nPairsToSkip
fractionAlignmentsPerPairInM = options.pairFrequency


# echo option choices/defaults
if not quiet:
   print "consensus_alignment  Sjölander lab, UC Berkeley"
   print " execution begun %s" % time.asctime()
   print "   Input seed pair:", seedpair
   print "   Input files:", files
   print "   Metric =", metric
   print "   Clustering method =", clusterMethod
   print "   Max clusters =", maxClusters
   print "   Distance threshold =", threshold
   if clusterMethod in ["heuristicHash", "dualPhase"]:
      print "   Hash on pairs of frequency approx. =", fractionAlignmentsPerPairInM
      print "   Hash on every", nPairsToSkip, "pair"
   print "   Consensus method =", consensusMethod
   if verbose:
      print "   Verbose mode"
   if quiet:
      print "   Quiet mode"
   print "   Reference alignment =", ref_alignment_file

#======================
# Read input alignments
#======================

# read reference alignment, if any

if ref_alignment_file:
   ref_alignment = Alignment()
   ref_alignment.readFromFile(ref_alignment_file)
else:
   ref_alignment = None

# read alignments to be clustered
# if filename specification includes *, ?, etc., must be in quotes
# so it's expanded here instead of by the shell.
if files:
  filelist = glob.glob(files)
elif seedpair:
  try:
    (seedX, seedY) = seedpair.split("_")
  except:
    raise ClusteringAndConsensusError, \
       "--seedpair value must be of form seedX_seedY"
filelist = matchmaker_seed_alignment_filenames(seedX, seedY, use_nr)
if len(filelist) == 0:
   raise ClusteringAndConsensusError,  "No input alignments."
if verbose:
   print "%d alignment files" % len(filelist)
   print "Reading alignments; printing one dot per 100:"
alignmentList = []
prev = None
_minPairs = 10     # default if not set to something else later
nalignments=0
discarded_alignments=0
for fname in filelist:
   a = Alignment()
   a.readFromFile(fname)
   if ref_alignment: a.reference = ref_alignment
   if prev != None:
      _minPairs = a.pairsNeededForCoverage(_coverageRequirement)
   if prev != None and not a.hasSameSeqsAs(prev):
      print >> sys.stderr, "Alignment %s omitted: sequences differ \
             from those in %s" % (a.fname, prev.fname)
   elif len(a.pairsNoGaps) <= _minPairs:
      # this alignment hardly aligns anything; not worth considering
      discarded_alignments = discarded_alignments + 1
      continue
   else:
      alignmentList.append(a)
      prev = a
   if verbose:
     nalignments = nalignments+1
     if (nalignments % 100) == 0:
       sys.stdout.write(".")
       sys.stdout.flush()
     if (nalignments % 5000) == 0:
       sys.stdout.write("%d", nalignments)
       sys.stdout.flush()

if verbose: sys.stdout.write("Done.\n")
     
if verbose and discarded_alignments:
   print "%d alignments with fewer than %d aligned pairs discarded" % \
      (discarded_alignments, _minPairs)

# check consistency of reference alignment with input alignments
if ref_alignment:
   if prev != None and not ref_alignment.hasSameSeqsAs(prev):
      print >> sys.stderr, "Sequence(s) in reference alignment differ(s) \
          from those in input alignment %s" % (ref_alignment.fname,
          prev.fname)

# get some info from first alignment in list
seed1 = alignmentList[0].name[0]
seed2 = alignmentList[0].name[1]

# Optional: print the distribution of distances among all pairs.
# testMetrics(alignmentList)

if verbose:
   print "First [up to] 10 input alignments:"
   i = 0
   for a in alignmentList:
      a.printLines()
      i = i + 1
      if i >=10: 
         break

#==================================
# Cluster by calling formClusters()
#==================================

clusters = formClusters(alignmentList, maxClusters, threshold,
      metric, clusterMethod, consensusMethod,
      fractionAlignmentsPerPairInM, nPairsToSkip,
      verbose)

#======================================
# Create and print consensus alignments
#======================================

print "Final number of clusters: %d" % len(clusters)
if consensusMethod != "none":
   alist_filename = alignment_consensus_csv_file(seed1, seed2)
   #alist_filename = seed1 + "_" + seed2 + "_consensus_alignments.csv"
   alist_file = open(alist_filename, "w")
   csv_filename = alignment_consensus_data_csv_file(seed1, seed2)
   #csv_filename = seed1 + "_" + seed2 + "_consensus_alignments_data.csv"
   csv_file = open(csv_filename, "w")
   print >> alist_file, "Alignment"
   print >> csv_file, csv_header()
  
for cluster in clusters:
   # compute consensus sequences if necessary
   # (with centroidJoining, cluster has already been computed and stored)
   if (clusterMethod != "centroidJoining") and (consensusMethod != "none"):
      clusterAlignmentSet = set()
      for i in cluster.members:
         clusterAlignmentSet.add(cluster.masterAlignmentList[i])
      cluster.consensus = consensusAlignment(clusterAlignmentSet,
              "minMedian", verbose)
   if ref_alignment:
      ref_distance = distanceBetweenAlignments(ref_alignment,
         cluster.consensus, metric)
      print "Distance from reference alignment: ", ref_distance
   if consensusMethod != "none":
      print >> csv_file, cluster._csv_line()
      print >> alist_file, os.path.join(cluster.consensus.dirname,
          cluster.consensus.fname)

tallyClusterSizes(clusters)

if ref_alignment:
   print "Reference alignment:"
   ref_alignment.printLines()
print "=============="

if not quiet:
   print " execution complete %s" % time.asctime()
