#!/bin/env python

"""
Input: matchmaker and certain nonmatchmaker alignment files
   matchmaker alignment files are ranked in a .csv file

Script builds a model for each alignment
  Models are built within the directory for each alignment.

TODO:
 Make this use build_one_model.py or build_models.py, under
    the control of a command-line option
 Command-line option for fast models
 Command-line option for which CSV file to get alignments from
"""

ohana = 0
salilab = 1

import os, glob, sys, re
if salilab:
  sys.path = ['/netapp/sali/ruchira/pylib/'] + sys.path
from matchmaker.shmm_shmm_lib import *
from matchmaker.build_model import *
from optparse import OptionParser

def compute_TSVMod_score(alignment_filename, pair,  score_csv_file):

    (seedX, seedY) = pair.split("_")
    fulldirname = os.path.dirname(alignment_filename)
    alignment_filename = os.path.basename(alignment_filename)

    basename = os.path.splitext(alignment_filename)[0]
    initial_dir = os.getcwd()
    os.chdir(fulldirname)
    pdb_filepath = os.path.join(single_seqs_dir(), seedX, "%s.pdb" % seedX[0:4])
    # link the template PDB file to this dir
    if os.path.exists("%s.pdb" % seedX[0:4]):
      os.system("rm %s.pdb" % seedX[0:4])
    os.system("ln -s %s ." % pdb_filepath)
    # run TSVMod
    os.system("perl /clusterfs/ohana/software/bin/TSVMod/score_model.pl " \
              + "-model %s.B99990001.pdb " % (seedY) \
              + "-alignment %s.ali -output_file %s.pred" % (basename, basename))
    # write scores to master list
    score_filename = "%s.pred_rmsd" % seedY
    if os.path.exists(score_filename):
      score_file = open("%s.pred_rmsd" % seedY, "r")
      score = float(score_file.readline().strip())
      score_file.close()
      score_csv_file.write("%s, %f\n" % (alignment_filename, score))
    else:
      score_csv_file.write("%s, None\n" % (alignment_filename))
    
    os.chdir(initial_dir)



def build(alignment_filename, pair, build=True, score=True, superpose=True,
          pdb_id='', job_prefix='model_'):

    # we have full pathname for nonmatchmaker alignment
    fulldirname = os.path.dirname(alignment_filename)
    alignment_filename = os.path.basename(alignment_filename)
    alignment_base_filename = os.path.splitext(alignment_filename)[0]
    model_base_filename = job_prefix + os.path.splitext(alignment_filename)[0]
    (seedX, seedY) = pair.split("_")
    pdb_filepath = os.path.join(single_seqs_dir(), seedX, "%s.pdb" % seedX[0:4])

    initial_dir = os.getcwd()
    os.chdir(fulldirname)
    # link the template PDB file to this dir
    os.system("rm -f %s.pdb" % seedX[0:4])
    os.system("ln -s %s ." % pdb_filepath)

    # create script to build model
    f = open(os.path.join(fulldirname, "%s.sh" % model_base_filename), "w")
    f.write("#!/bin/bash\n")
    # tell which node we're on?
    f.write("echo $HOSTNAME\n")
    # change into the alignment directory
    f.write("cd %s\n" % fulldirname)
    # build one standard (not fast) model
    if build:
      if ohana:
        f.write("build_one_model.py %s %s\n" % (fulldirname, pair))
      else:
        f.write("/netapp/sali/ModPipe/SVN/ext/mod/bin/modpy.sh ")
        f.write("python /netapp/sali/ruchira/bin/build_one_model.py %s %s\n" 
                % (fulldirname, pair))
    # run TSVMod
    if score:
      if ohana:
        f.write("perl /clusterfs/ohana/software/bin/TSVMod/score_model.pl ")
      else:
        f.write("perl /netapp/sali/ruchira/TSVMod/main/score_model.pl ")
      f.write("-model %s.B99990001.pdb " % (seedY) \
                + "-alignment %s.ali " % alignment_base_filename \
                + " -output_file %s.pred\n" % alignment_base_filename)

    if superpose:
      if salilab:
        f.write("/netapp/sali/ModPipe/SVN/ext/mod/bin/modpy.sh python ")
        f.write("/netapp/sali/ruchira/bin/")
      f.write("superpose_model_on_true_structure.py ")
    f.write("%s %s %s.B99990001.pdb %s\n" \
            % (pdb_id, seedY, seedY, alignment_filename))

    f.write("echo '%s is done' > %s/done\n" % (alignment_filename, fulldirname))
    f.close()

    # submit job to queue
    if ohana:
      os.system("qsub -j oe -o %s.out %s.sh" % (model_base_filename, model_base_filename))
    else:
      os.system("qsub -S /bin/bash -cwd -j y -o %s.out %s.sh >& qsub.out" % (model_base_filename, model_base_filename))
    os.chdir(initial_dir)

def main():
  parser = OptionParser()
  parser.add_option("--build", dest="build", action="store_true", default=True,
      help="Build models (default is to build)")
  parser.add_option("--nobuild", dest="build", action="store_false",
      default=True,
      help="Do not build models (default is to build)")
  parser.add_option("--score", dest="score", action="store_true", default=True,
      help="Score models with TSVMod (default is to score)")
  parser.add_option("--noscore", dest="score", action="store_false",
      default=True,
      help="Do not score models with TSVMod (default is to score)")
  parser.add_option("--superpose", dest="superpose", action="store_true", 
      default=True,
      help="Superpose models on true structure (default is to superpose)")
  parser.add_option("--nosuperpose", dest="superpose", action="store_false",
      default=True,
      help="Do not uperpose models on true structure (default is to superpose)")
  parser.add_option("--job_prefix", dest="job_prefix", default="model_",
      help='Prefix for the job scripts (default: "model_")')
  parser.add_option("--pdb_id", dest="pdb_id", default="",
      help="PDB ID of the true structure")
  (options, args) = parser.parse_args()
  if len(args) < 1:
    print "Usage: %s [options] <pair>" % sys.argv[0]
    sys.exit(0)

  pair = args[0]
  (seedX, seedY) = pair.split("_")
  pdb_id = options.pdb_id
  if options.superpose and pdb_id == '':
    pdb_id = seedY
  # can't write scores to file because they won't be ready until after
  # queued jobs are done
  #score_csv_filename = os.path.join(align_dir(seedX, seedY), "tsvmod_rmsd.csv")
  #score_csv_file = open(score_csv_filename, "w")

  for old_alignment_filename in \
         matchmaker_seed_alignment_filenames(seedX, seedY, use_nr=True):

    alignment_filename = \
        old_alignment_filename.replace('/clusterfs/ohana/bpg/matchmaker',
                                    shmm_root_dir())
    build(alignment_filename, pair, build=options.build, score=options.score,
          superpose=options.superpose, pdb_id=pdb_id, job_prefix=options.job_prefix)
    #compute_TSVMod_score(alignment_filename, pair,  score_csv_file)

  #score_csv_file.close()

  # sort results
  # cmd = 'sort -r -n -t "," --key=2 %s > %s.sorted' % (score_csv_filename,
      #score_csv_filename)
  #os.system(cmd)

if __name__ == "__main__":
  main()
