#!/bin/env python

# Within MODELLER, build models.  Argument is directory containing one alignment
# file (assumed to have extension "afa")

ohana = 1
salilab = 0

import csv, glob, os, re, sys
if salilab:
  sys.path = ['/netapp/sali/ruchira/pylib/'] + sys.path
from modeller import *
import matchmaker.build_model

# First argument is alignment directory.
if len( sys.argv ) == 1 :
   print "Need directory of alignment files on command line"
   sys.exit( 1 )

alignment_dir = sys.argv[1]
if not os.access( alignment_dir, os.F_OK ) :
   print "Could not access", alignment_dir
   sys.exit( 1 )

# Second argument is template/target pair
pair = sys.argv[2]

# Third, optional argument says to do a very fast model
fast_model = (len( sys.argv ) > 2)

# Get alignment file.
alignments = glob.glob( "%s/*.afa" % alignment_dir )
if not alignments :
   print >> sys.stderr, "Did not find alignment files (.afa) in directory %s" % \
       alignment_dir
   sys.exit( 1 )
if len(alignments) > 1:
  print >> sys.stderr, "More than one alignment file found in %s; using first" % \
       alignment_dir

# Build model from
# the alignment.  Write to results file in directory.

results_file = "%s/model_quality_results.csv" % alignment_dir
fp_results = open( results_file, "w" )
fp_results.write("Alignment,DOPE,GA341\n")
alignment = alignments[0]
dope_score, ga341_score = matchmaker.build_model.build_model( alignment, pair,
    fast_model )
fp_results.write( "%s,%s,%s\n" % ( alignment, dope_score, ga341_score ) )

# Done.
fp_results.close()

