#!/bin/env python

"""
Get PDB and sequence files for a template.
"""

import os, sys, glob
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <template_id>" % sys.argv[0]
  sys.exit(1)
template = sys.argv[1]

# go into the appropriate subdir of single_seqs
dir = os.path.join(single_seqs_dir(), template)
make_dir_exist(dir)
initial_dir = os.getcwd()
os.chdir(dir)

# get pdb file
cmd = "wget www.rcsb.org/pdb/files/%s.pdb" % template[0:4]
print cmd
os.system(cmd)

# get sequence file and change header
cmd = "fastacmd -d $BLASTDB/pdb_rcsb -s %s_%s | sed '/>/c\>%s'> %s.fa" % \
    (template[0:4], template[4], template, template)
print cmd
os.system(cmd)

# go back to the  dir we started from
os.chdir(initial_dir)
