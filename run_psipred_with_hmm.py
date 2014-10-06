#!/usr/bin/python

"""
Input: A profile HMM in SAM (.mod) format
Output: PSIPRED secondary structure prediction
"""

import os, sys
from optparse import OptionParser

psipred_data_dir = "/clusterfs/ohana/software/lib/psipred/data"

if len(sys.argv) < 3:
  print "Usage: %s runname hmm_filepath" % sys.argv[0]
  sys.exit(0)

runname = sys.argv[1]
hmm_filepath = sys.argv[2]

cmd = "sam2psi %s -i %s" % (runname, hmm_filepath)
print cmd
os.system(cmd)

# we must run makemat on a copy of the .ckp file, because
# it will overwrite the original
cmd = "cp %s.ckp %s.makemat.ckp" % (runname, runname)
print cmd
os.system(cmd)

cmd = "echo %s.makemat.ckp > %s.pn" % (runname, runname)
print cmd
os.system(cmd)

cmd = "echo %s.cks > %s.sn" % (runname, runname)
print cmd
os.system(cmd)

cmd = "makemat -P %s" % (runname)
print cmd
os.system(cmd)

# the name of the makemat output file is stored in a file
makemat_matrix_record_filename = runname + ".mn"
makemat_matrix_record_file = open(makemat_matrix_record_filename, "r")
makemat_matrix_filename = makemat_matrix_record_file.readline().strip()
print makemat_matrix_filename
cmd = "cp %s %s.mtx" % (makemat_matrix_filename, runname)
print cmd
os.system(cmd)


cmd ="psipred %s.mtx %s/weights.dat %s/weights.dat2 %s/weights.dat3 %s/weights.dat4 > %s.ss" % \
   (runname, psipred_data_dir, psipred_data_dir, psipred_data_dir,
    psipred_data_dir, runname)
print cmd
os.system(cmd)

cmd = "psipass2 %s/weights_p2.dat 1 1.0 1.0 %s.ss2 %s.ss > %s.horiz" % \
    (psipred_data_dir, runname, runname, runname)
print cmd
os.system(cmd)
