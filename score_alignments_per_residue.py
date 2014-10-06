#!/usr/bin/env python

from optparse import OptionParser
import os, sys, glob, string
from matchmaker.alignments_lib import *
from matchmaker.cluster_consensus_lib import _GAP
from matchmaker.shmm_shmm_lib import *

#====================
# Process command line
#====================
parser = OptionParser()

"""
parser.add_option("-f", "--files", dest="files",
      help="files containing input alignments; enclose in quotes if using *, ?, etc.", metavar="FILES")
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
if not options.files:
   parser.error("Option -f required")
files = options.files
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
   print "consensus_alignment  Sjolander lab, UC Berkeley"
   print " execution begun %s" % time.asctime()
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

"""



if len(sys.argv) < 2:
  print "Usage: %s <pair> [<N>]" % sys.argv[0]
  sys.exit(0)

pair = sys.argv[1]
if len(sys.argv) > 2:
  max_alignments = int(sys.argv[2])
else:
  max_alignments = 5

(alignment_pathnames, model_pathnames, alignment_csv_lines) = \
  top_models_by_pred_rmsd(pair, max_alignments)
num_alignments = len(alignment_pathnames)

selected_alignment_filenames = [os.path.basename(file) for file in
    alignment_pathnames]

class ColumnInfo:
  """ keeps track of best score, total score in a column of CSV file"""
  def __init__(self):
    self.best_score = -1.0
    self.best_i = -1
    self.total_score = float(0)
    self.scoreDict = {}

columnInfoDict = {}
alignmentDict = {}
maxj = -1

# for each alignment, tally the aligned pairs
for i in range(0,num_alignments):
  alignment_filename = alignment_pathnames[i]
  alignment_basename = selected_alignment_filenames[i]
  # get the pred_no35 (don't bother with  other scores)
  line = alignment_csv_lines[i]
  fields = line.rstrip().split(',')
  pred_no_35 = string.atof(fields[2])
  # read the alignment into an Alignment object
  alignment = Alignment()
  alignmentDict[alignment_basename] = alignment
  alignment.readFromFile(alignment_filename)
  alignment.positionPairsFromLines()
  # tally the aligned pairs in this alignment
  for (i, j) in alignment.pairsInclGaps:
    if j != _GAP:
      if j > maxj:
        maxj = j
      if j not in columnInfoDict or columnInfoDict[j] == None:
        columnInfoDict[j] = ColumnInfo()
      if i in columnInfoDict[j].scoreDict and columnInfoDict[j].scoreDict[i] != None:
        columnInfoDict[j].scoreDict[i] = columnInfoDict[j].scoreDict[i] + pred_no_35
      else:
        columnInfoDict[j].scoreDict[i] = pred_no_35
      columnInfoDict[j].total_score = columnInfoDict[j].total_score + pred_no_35
      if columnInfoDict[j].scoreDict[i] > columnInfoDict[j].best_score:
        columnInfoDict[j].best_i = i
        columnInfoDict[j].best_score = columnInfoDict[j].scoreDict[i]

# get aligned sequences from the first alignment
targetseq = alignmentDict[selected_alignment_filenames[0]].seq[1]
templateseq = alignmentDict[selected_alignment_filenames[0]].seq[0]

# write first part of header line
sys.stdout.write('Tgt Res Num, Tgt Res,')
sys.stdout.write(' Most Likely Tmplt Res Num, Most Likely Tmplt Res, Score Most Likely Tmplt Res')
# for all alignments named on command line, add  fields to header
for alignment_filename in selected_alignment_filenames:
  aerated_name = alignment_filename.replace("_"," ").replace("."," ")
  sys.stdout.write(',"Tmplt Res Num In %s","Tmplt Res Num"' % aerated_name)
  sys.stdout.write(',"Score"')
sys.stdout.write("\n")

# for each alignment position, write a row to the table
for j in range(0, maxj + 1):
  sys.stdout.write("%d,%s" % (j, targetseq[j]))
  best_score = columnInfoDict[j].best_score / columnInfoDict[j].total_score
  if columnInfoDict[j].best_i == _GAP:
    sys.stdout.write(",-1,-,%0.3f" % best_score)
  else:
    sys.stdout.write(",%d,%s,%0.3f" % (columnInfoDict[j].best_i,
                                  templateseq[columnInfoDict[j].best_i],
                                  best_score))
  for alignment_filename in selected_alignment_filenames:
    alignment = alignmentDict[alignment_filename]
    i = alignment.seqpairs[1][j]
    score = columnInfoDict[j].scoreDict[i] / columnInfoDict[j].total_score
    if i == _GAP:
      sys.stdout.write(",-1,-,%0.3f" % score)
    else:
      sys.stdout.write(",%d,%s,%0.3f" % (i, templateseq[i], score))
  sys.stdout.write("\n")
