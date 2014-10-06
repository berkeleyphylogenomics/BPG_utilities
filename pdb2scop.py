#!/usr/bin/python
import sys
from BPG_common.fasta import *

"""
May 3, 2008   Terry Farrah   Sjolander Lab, UC Berkeley
Retrieves the SCOP domain identifier from the Astral PDB40 database
for a given sequence ID.
"""


def pdb2scop(seqid):
  # use fastacmd to retrieve the sequence record from the database,
  # putting the record in a temporary file
  seqfilename = seqid + ".fa"
  pdb40filename = "/home/ruchira/CASP8/astral40.fasta"
  cmd = "fastacmd -s %s -d %s >| %s" % (seqid, pdb40filename, seqfilename)
  os.system(cmd)

  # read the file using a function from the local fasta library 
  seq_list = ReadSequencesList(seqfilename)
  if len(seq_list) == 0:
    print >> sys.stderr, "seqID %s not found in %s" % (seqid, pdb40filename)
    sys.exit(0)
  seq_tuple = seq_list[0]

  # get the scop domain from the  header
  header = seq_tuple[0]
  #print header
  domain = header.split(" ")[1]

  # clean up: remove sequence file
  os.system("rm %s" % (seqfilename))

  # return the domain as a string
  return(domain)

# Command line interface
# get sequence ID from command line
def main():
  if len(sys.argv) < 2:
    print "Usage: pdb2scop <pdbid>"
    sys.exit(0)
  seqid = sys.argv[1]

  # call pdb2scop()
  pdb = pdb2scop(seqid)

  # print the result
  print pdb

if __name__ == '__main__':
     main()
