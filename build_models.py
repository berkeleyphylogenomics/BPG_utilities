#!/usr/bin/python
# Within MODELLER, build models.  Argument is directory of alignment files
# (assumed to have extension "afa") and csv file of alignment files and CS 
# scores (assumed to have extension "csv"), which may include additional 
# alignments).

# 4/4/08: large sections commented out so that (a) script assumes
# alignments are already in alignment dir, and (b) script does not
# look for or output CS scores. Did this because I modified the
# previous script to already link the alignments to the alignment
# dir, and  also disabled CSV file writing due to a bug I didn't find.
# 6/3/08: restored those sections.

import csv, glob, os, re, sys
from matchmaker.shmm_shmm_lib import *
from modeller import *
import build_model

debug1 = 0
debug2 = 0  # Stop after this many (0 = don't stop)

# First argument is alignment directory.

if len( sys.argv ) < 2:
   print "Usage: %s <pair_id>" % (sys.argv[0])
   sys.exit( 1 )

pair_id = sys.argv[1]
seeds = pair_id.split("_")
if len(seeds) != 2:
  print >> sys.stderr, "%s: pair %s not of form seedX_seedY" % \
     (sys.argv[0], pair_id)
seedX = seeds[0]
seedY = seeds[1]
alignment_dir = selected_align_dir(seedX, seedY)
if not os.path.exists(alignment_dir):
  print >> sys.stderr, "%s: alignment directory %s does not exist" % \
      (sys.argv[0], alignment_dir)
  sys.exit(1)

# Second, optional argument says to do a very fast model

fast_model = (len( sys.argv ) > 1)


# Look for csv file of alignments and CS scores.  This will also pick up the
# results file, so discard those.

csv_files = glob.glob( "%s/*.csv" % alignment_dir )
alignment_cs_files = []
for csv_file in csv_files :
   if not re.search( "model", csv_file ) :
      alignment_cs_files.append( csv_file )

if not alignment_cs_files :
   print "Could not find csv file of alignments,CS scores"
   sys.exit( 1 )

n_files = len( alignment_cs_files )
if n_files != 1 :
   print "Found %s csv files.  Using %s" % ( n_files, alignment_cs_files[0] )

alignment_cs_file = alignment_cs_files[0]

# Read csv file.  Save CS scores as dictionary array indexed by alignment file
# name.

alignment_cs_fp = open( alignment_cs_file, "r" )
reader = csv.reader( alignment_cs_fp )

# if any/all CS score(s) do/does not exist, value is None
cs = {}
# discard header row
header_row = reader.next()
for row in reader:
   cluster_size = int(row[0])
   file = row[1]
   cs_i = row[2]
   alignment_name = file.replace( "_CS.out", "" )
   cs[alignment_name] = cs_i


# Get list of alignments.


alignments = glob.glob( "%s/*.afa" % alignment_dir )
if not alignments :
   print "Did not find alignment files (.afa) in directory %s" % alignment_dir
   sys.exit( 1 )

# Check that each alignment has a CS score.

for alignment in alignments :
   dir, alignment = os.path.split( alignment )
   alignment_name = alignment.replace( "_original.afa", "" )
   if alignment_name not in cs.keys() :
      print "Could not find CS score for alignment %s in file %s" \
                                              % ( alignment, alignment_cs_file )

# Build two models (seedX as template with seedY as target, and vice-versa) from
# each alignment.  Write, along with CS score, to results file in directory.

results_file = "%s/model_quality_results.csv" % alignment_dir
fp_results = open( results_file, "w" )
i_alignment = 0
for alignment in alignments :
   dope_score_seedX_template, ga341_score_seedX_template, \
   dope_score_seedY_template, ga341_score_seedY_template \
             = build_model.build_model2( alignment, fast_model )

   # Get CS score.

   dir, alignment = os.path.split( alignment )
   alignment_name = alignment.replace( "_original.afa", "" )
   cs_score = cs.get( alignment_name, -99999.9 )
   fp_results.write( "%s,%s,%s,%s,%s,%s\n" % ( alignment, \
                                               cs_score, \
                                               dope_score_seedX_template, \
                                               ga341_score_seedX_template, \
                                               dope_score_seedY_template, \
                                               ga341_score_seedY_template ) )
   i_alignment += 1
   if debug2 :
      if i_alignment >= debug2 :
         break

# Done.

fp_results.close()

