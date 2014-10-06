#!/usr/bin/python

"""
Converts a DSSP secondary-structure description 
to a three-state description in psipred format..

Uses the following correspondence:
  --------------------------
  | DSSP  |H|G|I|E|B|T|S|''|
  --------------------------
  |3-state|H|H|H|E|E|L|L| L|
  --------------------------
This table appears at
http://cubic.bioc.columbia.edu/eva/doc/measure_sec.html
"EVA measures for secondary structure protein accuracy"
Author: Burkhard Rost rost@columbia.edu
However, here 'C' is used in place of 'L', 
in accordance with the PSIPRED convention.
"""
import os, sys
import BPG_common.fasta
from optparse import OptionParser

dssp_to_threestate_dict = {}
dssp_to_threestate_dict['H'] = 'H'
dssp_to_threestate_dict['G'] = 'H'
dssp_to_threestate_dict['I'] = 'H'
dssp_to_threestate_dict['E'] = 'E'
dssp_to_threestate_dict['B'] = 'E'
dssp_to_threestate_dict['T'] = 'C'
dssp_to_threestate_dict['S'] = 'C'
dssp_to_threestate_dict[ ''] = 'C'

def dssp_to_threestate(dssp_pred):
  if dssp_pred in dssp_to_threestate_dict:
    threestate = dssp_to_threestate_dict[dssp_pred]
    if threestate == None:
      return ""
    else:
      return threestate
  return ""

def main():
  parser = OptionParser()
  parser.add_option("-i", "--input",
    dest="dssp_filename",
    metavar="FILE",
    help="dssp output file")
  parser.add_option("-s", "--seqfile",
    dest="seq_filename",
    metavar="FILE",
    help="sequence file in fasta format, with X's for residues lacking coords (optional)")
  parser.add_option("-o", "--output",
    dest="out_file",
    metavar="FILE",
    help="psipred vformat output file")
  (options, args) = parser.parse_args()
  # check that necessary options are given, and that values are valid
  # assign option values to variables
  if not options.dssp_filename:
    parser.error("Option -i required")
  dssp_filename = options.dssp_filename
  seq_filename =  options.seq_filename
  if not options.out_file:
    out_file = sys.stdout
  else:
    out_file = open(options.out_file, "w")

  f = open(dssp_filename)
  lines = f.readlines()
  f.close()
  if seq_filename:
    if seq_filename.endswith(".mask"):
      f = open(seq_filename, "r")
      seq = f.readline().strip()
      def is_not_star(c): return (c != "*")
      seqlen = len(filter(is_not_star, seq))
    else:
      seq = BPG_common.fasta.ReadOneSequence(seq_filename)
      seqlen = len(seq)

  out_file.write("# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n\n")
  found_description_start = False
  vector_of_threestate = {}
  vector_of_threestate['C'] = "  1.000  0.000  0.000"
  vector_of_threestate['H'] = "  0.000  1.000  0.000"
  vector_of_threestate['E'] = "  0.000  0.000  1.000"
  vector_of_threestate['X'] = "  0.333  0.333  0.333"

  # open DSSP file
  # open sequence file
  # there should be a one-to-one correspondence between
  # DSSP file records and non-X sequence characters
  line_num = 0
  res_num = 0
  do_pad = False
  # process all residue records in DSSP file
  for line in lines:
    line_num = line_num + 1
    line = line.rstrip()
    if line[2] == '#':
      found_description_start = True
      if seq_filename:
	if len(lines) - line_num == seqlen - 1:
	  print >> sys.stderr, \
	     "WARNING: DSSP descriptor has one fewer res than seq; will pad output at end"
	  do_pad = True
	elif len(lines) - line_num != seqlen:
	  print >> sys.stderr, \
	     "ERROR: Lengths of DSSP description (%d), sequence(%d), differ" \
	     % ( len(lines) - line_num, seqlen)
	  sys.exit(1)
    
    elif found_description_start:
      if len(line) >= 17:
        #num = line[0:5].strip()
        residue_num = line[5:10].strip()
        # For 1kl9A, Terry found a DSSP output line that
        # was missing a residue #. The line did not
        # have any corresponding ATOM records in the
	# PDB file, so it makes sense to ignore it.
        try:
          dssp_res_num = int(residue_num)
        except:
          continue
        
        res_num = res_num + 1
        if seq_filename:
          seq_char = seq[res_num-1]
          while seq_char == '*':
            out_file.write("%4s X X %s\n" % (res_num, vector_of_threestate['X']))
            res_num = res_num+1
            seq_char = seq[res_num-1]

        aa = line[12:14].strip()
        if seq_filename:
          if aa.upper() != seq[res_num-1].upper():
	    print >> sys.stderr, \
	       "WARNING: Residue at pos %d in dssp (%s) doesn't match in seq (%s)" \
	       % ( res_num, aa, seq[res_num-1])
 
        dssp_pred = line[14:17].strip()
        threestate = dssp_to_threestate(dssp_pred)
        if threestate != "":
          out_file.write("%4d %s %s %s\n" % (res_num, aa, threestate,
                                      vector_of_threestate[threestate]))

  # write records for any asterisks remaining in sequence.
  if seq_filename:
    res_num = res_num + 1
    while res_num <= len(seq):
      seq_char = seq[res_num-1]
      if seq_char != '*':
        print >> sys.stderr, \
           "WARNING: Residue at pos %d in seq (%s) has no dssp record; X written." \
           % ( res_num, seq_char )
      out_file.write("%4s X X %s\n" % (res_num, vector_of_threestate['X']))
      res_num = res_num+1

  # if sequence has more chars than DSSP file has lines, pad with an  X
  if do_pad:
    res_num = res_num + 1
    out_file.write("%4s %s %s %s\n" % (res_num, 'X', 'X',
                                      vector_of_threestate['X']))
  out_file.close()

  
if __name__ == '__main__':
  main()
