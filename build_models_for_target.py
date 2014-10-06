#!/usr/bin/python

import os, glob, sys, re

if len(sys.argv) < 1:
  print "Usage: %s <pair>" % sys.argv[0]
  sys.exit(0)

# 05/20/2008 TODO: Make this work for alignments other than SHMM-SHMM
# profile-profile using YL scoring, under the control of a command-line option
# 05/20/2008 TODO: Make this use build_one_model.py or build_models.py, under
# the control of a command-line option
basedir = "/home/ruchira/SHMM-SHMM/current_method/results"
def build(alignment):
    fields = alignment.split('_')
    seedX = fields[0]
    seedY = fields[1]
    hmmX = fields[2]
    hmmY = fields[3]
    fulldirname = os.path.join(basedir, "%s_%s_align/%s/%s/YL" 
                                          % (seedX, seedY,
                                              hmmX, hmmY))
    os.system("cp /home/ruchira/src/research/modelling/build_model.py %s" \
      % fulldirname)
    os.system("cp /home/ruchira/src/research/modelling/build_one_model.py %s" \
      % fulldirname)
    modelbasename = "%s_%s_%s_%s_YL_model" % (seedX, seedY, hmmX, hmmY)
    f = open(os.path.join(fulldirname, "%s.sh" % modelbasename), "w")
    f.write("#!/bin/bash\n")
    f.write("echo $HOSTNAME\n")
    f.write("cd %s\n" % fulldirname)
    f.write("python %s/build_one_model.py %s\n" \
      % (fulldirname, fulldirname))
    f.write("echo '%s is done' > %s/done\n" % (alignment, fulldirname))
    f.close()
    cur_dir = os.getcwd()
    os.chdir(fulldirname)
    os.system("qsub -j oe -o %s.out %s.sh" % (modelbasename, modelbasename))
    os.chdir(cur_dir)

csvs = glob.glob("%s/%s_align/*top*.csv" % (basedir, sys.argv[1]))
for csv in csvs:
  f = open(csv)
  f.readline()
  lines = f.readlines()
  f.close()
  for line in lines:
    fields = line.split(',')
    build(fields[0])
