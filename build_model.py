# Within MODELLER environment, setup to build homology model, etc.

import os, sys, tempfile, string, re

from modeller import *
from modeller.automodel import *
from shmm_shmm_lib import *

# 06/12/08 TODO: Make removal of the initial methionine a command-line option

debug1 = 0

bpg = 0
ohana = 0
salilab = 1


if bpg :
   pdb_dir = "/home/bpg/pdb_files"
   pdb_prefix = ""
   pdb_ext = "pdb"
elif ohana:
   pdb_dir = "/clusterfs/ohana/external/pdb_SHMM-SHMM"
   pdb_prefix = "pdb"
   pdb_ext = "ent"
else :
   # Salilab 
   pdb_dir = "/netapp/database/pdb/remediated/pdb"
   pdb_prefix = "pdb"
   pdb_ext = "ent.gz"


# ------------------------------------------------------------------------------

def find_best_score( models, score_name ) :
   if score_name == "GA341 score" :
      fac = 1.0
   else :
      fac = -1.0

   if models :
      best_pos_score = -1.0e8
      for model in models :
         if score_name == "GA341 score" :
            score = model[score_name][0]
         else :
            score = model[score_name]

         if debug1 :
            print "[find_best_score] score:", score

         score = float( score )
         pos_score = fac * score
         if pos_score > best_pos_score :
            best_pos_score = pos_score
            best_model = model["name"]

      best_score = fac * best_pos_score
         
   else :

      # Missing data.

      best_score = -99999.9
      best_model = "No_successful_model_build"

   return best_model, best_score


# ------------------------------------------------------------------------------
def pdb_filename( pdb_id ) :
   if salilab:
     return "%s/%s/%s%s.%s" % ( pdb_dir, pdb_id[1:3], pdb_prefix, pdb_id, pdb_ext )
   else:
     return "%s/%s%s.%s" % ( pdb_dir, pdb_prefix, pdb_id, pdb_ext )


# ------------------------------------------------------------------------------
def build_model2( seedX_seedY_alignment_file, fast_model=False ) :
   # 7/16/08: still needs to be adapted to compressed alignment files

   # Build two models from alignment: one with seedX as template and seedY as
   # target, and vice-versa.

   # First: just use file as-is.

   dope_score_seedX_template, ga341_score_seedX_template \
        = build_model(seedX_seedY_alignment_file, None, fast_model=fast_model)

   # Second: create a temporary alignment file, reversing the order.

   lines = open( seedX_seedY_alignment_file, "r" ).readlines()

   dir, alignment_file = os.path.split( seedX_seedY_alignment_file )
   basename, ext = os.path.splitext( alignment_file )
   seedY_seedX_alignment_file = tempfile.mktemp( basename )
   fp = open( seedY_seedX_alignment_file, "w" )
   fp.write( lines[2] )
   fp.write( lines[3] )
   fp.write( lines[0] )
   fp.write( lines[1] )

   fp.flush()
   fp.close()

   dope_score_seedY_template, ga341_score_seedY_template \
                                     = build_model(seedY_seedX_alignment_file,
                                                    None, fast_model=fast_model)

   return dope_score_seedX_template, ga341_score_seedX_template, \
          dope_score_seedY_template, ga341_score_seedY_template \

methionine_re = re.compile("(-*)M")
# ------------------------------------------------------------------------------
def build_model( template_target_afa_file, pair, fast_model=False ) :

   # Argument: template-target alignment file name.
   # Aligned FASTA format.  We'll assume sequences are on one line, and that 
   # headers look like ">1cc8A".
   # Template first, then target.

   if not os.access( template_target_afa_file, os.F_OK ) :
      print "Could not access %s" % template_target_afa_file
      return

   # Read file.

   """
   # 7/15/08: replace readlines() with function that can handle compressed files
   afa_lines = open( template_target_afa_file ).readlines()

   # Get template and target IDs, without ">" or EOL.
   # >1cc8a
   # 012345

   template_pdb_id_chain_id = afa_lines[0].strip()
   template_pdb_id_chain_id = template_pdb_id_chain_id[1:]

   target_id = afa_lines[2].strip()
   target_id = target_id[1:]

   """

   (template_pdb_id_chain_id, template_aligned_seq), (target_id,
     target_aligned_seq) = read_alignment_file(template_target_afa_file)
   template_pdb_id = template_pdb_id_chain_id[0:4]
   template_chain_id = template_pdb_id_chain_id[4:5]  # Works, even if length 4.
   if debug1 :
      print "template_pdb_id %s, template_chain_id %s" \
                                        % ( template_pdb_id, template_chain_id )
   if debug1 :
      print "target_id:", target_id
   

   # Create PIR-format file.

   dir, alignment_file = os.path.split( template_target_afa_file )
   basename, ext = os.path.splitext( alignment_file )
   tmpname = tempfile.mktemp( basename )
   template_target_pir_file = basename + ".ali"
   if debug1 :
      print "template_target_pir_file:", template_target_pir_file

   fp = open( template_target_pir_file, "w" )
   fp.write( ">P1;%s\n" % ( template_pdb_id_chain_id ) )
   structure_line \
         = "structureX:%s:FIRST:%s:LAST:%s:undefined:undefined:-1.00:-1.00\n" \
                     % ( template_pdb_id, template_chain_id, template_chain_id )
   fp.write( structure_line )
   #sequence_line = afa_lines[1].strip() + "*\n"
   sequence_line = template_aligned_seq + "*\n"
   sequence_line = sequence_line.upper()
   sequence_line = string.replace(sequence_line, '.', '-')
#   methionine_match = methionine_re.match(sequence_line)
#   if methionine_match != None:
#    sequence_line = string.replace(sequence_line, 'M', '-', 1)
   fp.write( sequence_line )

   fp.write( ">P1;%s\n" % target_id )
   fp.write( "sequence:%s:    : :     : ::: 0.00: 0.00\n" % target_id )
   #sequence_line = afa_lines[3].strip() + "*\n"
   sequence_line = target_aligned_seq + "*\n"
   sequence_line = sequence_line.upper()
   sequence_line = string.replace(sequence_line, '.', '-')
#   methionine_match = methionine_re.match(sequence_line)
#   if methionine_match != None:
#    sequence_line = string.replace(sequence_line, 'M', '-', 1)
   fp.write( sequence_line )

   fp.flush()
   fp.close()

   # Make sure template PDB file is available.

   template_pdb_file = pdb_filename( template_pdb_id )
   if not os.access( template_pdb_file, os.F_OK ) :

      # Not there.  If on bpg, run download script.

      if bpg :
         cmd = "source /etc/profile;\n"
         cmd += "/var/www/html/book/get_pdb_file.php %s " % template_pdb_id
         if debug1 :
            print cmd

         os.system( cmd )

   # Error if still not there.

   if not os.access( template_pdb_file, os.F_OK ) :
      print "Could not access %s" % template_pdb_file
      sys.exit(0)

   # ...........................................................................
   # Instance of MODELLER environment.

   env = environ()

   env.io.atom_files_directory = pdb_dir

   knowns = template_pdb_id + template_chain_id
   a = automodel( env, 
                  alnfile=template_target_pir_file,
                  knowns=knowns, 
                  sequence=target_id,
                  assess_methods=( assess.DOPE,
                                   assess.GA341 
                                 )
                )
   if fast_model:
      a.very_fast()
      a.starting_model = 1
      a.ending_model = 1
      #a.final_malign3d = True    # shown in example from manual
   else:
      a.starting_model = 1
      a.ending_model = 5

   a.make()

   # ...........................................................................
   # Get a list of all successfully built models from a.outputs

   ok_models = filter(lambda x: x['failure'] is None, a.outputs)

   # Pull out the best DOPE score (most negative) and best GA341 score
   # (highest).

   best_DOPE_model, best_DOPE_score = find_best_score( ok_models, "DOPE score" )
   best_GA341_model, best_GA341_score = find_best_score( ok_models, \
                                                                 "GA341 score" )

   # Return results.

   return best_DOPE_score, best_GA341_score
