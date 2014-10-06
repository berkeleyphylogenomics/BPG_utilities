#!/usr/bin/python

"""
May 7, 2008  Terry Farrah   Sjolander Lab, UC Berkeley

Input: DSSP output file (secondary structure for a PDB file)
Output: a fasta format sequence file, with a secondary
  structure code letter at each position, one for each
  amino acid
"""

import sys

def main():
  # get dssp file from command line and open
  if len (sys.argv) < 2:
    print >> sys.stderr, "Usage: %s <dssp-output-file>" % \
        sys.argv[0]
    sys.exit(0)
  dssp_filename = sys.argv[1]
  dssp_file = open(dssp_filename, "r")

  # add a .fa extension for the output filename and open
  out_filename = dssp_filename + ".fa"
  out_file = open(out_filename, "w")

  # put doc ine into output file
  print >> out_file, "> %s secondary structure string" % dssp_filename

  # read lines until find one with # in third position
  line = dssp_file.readline()
  while line[2] != "#":
    line = dssp_file.readline()

  # after that, read one residue per line. 17th char,
  charcount = 1
  while line != "":
    sschar = line[16]

    # if not space, is SS. If space, put X for undetermined.
    if sschar == " ":  sschar = "X"
    # write char to ouptut.
    out_file.write(sschar)
    if charcount % 60 == 0:
      out_file.write("\n")
    
    line = dssp_file.readline()
    charcount = charcount + 1

  # write newline and close files.
  out_file.write("\n")
  out_file.close()
  dssp_file.close()


if __name__ == "__main__":
  main()
