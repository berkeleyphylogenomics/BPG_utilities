#!/usr/bin/python

import os, sys

"""
May 6, 2008   Terry Farrah
Input: hmm and corresponding seed sequence. PDB db available.
Output: msa of significant PDB matches to seed

Initial implementation, May 6 2008:
This script must be run in the same directory as the hmm.
"""
pdb_db_filename = "$BLASTDB/pdb_rcsb"
fastacmd_in_filename = "fastacmd.in"
fastacmd_out_filename = "fastacmd.out"
eval_threshold = -2

def parse_eval(eval_string):
  '''
  Return a (coefficient, exponent) tuple for a text eval of form 2.1e-04
  '''
  if eval_string.find("e") > -1:
    eval = eval_string.split("e")
    coefficient = float(eval[0])
    exponent = int(eval[1])
  else:
    coefficient = float(eval_string)
    exponent = 0
  return coefficient, exponent


def get_significant_hits_from_hmmscore_output(hmmscore_out_filename, threshold_exponent):
  ''' 
  Return seq IDs for hits with eval exponent <= threshold_exponent.
  '''
  hit_list = []
  hmmscore_out_file = open(hmmscore_out_filename, "r")
  hmm_output_line = hmmscore_out_file.readline()
  # skip comment lines
  while hmm_output_line != "" and hmm_output_line[0] == "%":
    hmm_output_line = hmmscore_out_file.readline()
  while hmm_output_line != "":
    elements = hmm_output_line.split()
    [subject, length, simple, reverse, hmm_eval] = \
        hmm_output_line.split()[:5]
    hmmscore_eval_coefficient, hmmscore_eval_exponent = parse_eval(hmm_eval)
    if hmmscore_eval_exponent <= threshold_exponent:
      hit_list.append(subject)
    hmm_output_line = hmmscore_out_file.readline()
  return ( hit_list )


def main():

  if len(sys.argv) < 2:
    print >> sys.stderr, "Usage: %s <hmm_filename>" % (sys.argv[0])
    print >> sys.stderr, "  NOTE: hmm_filename must be of the form"
    print >> sys.stderr, "  seed.ext (e.g. seed.mod, seed.ghmm)"
    print >> sys.stderr, "  and a fasta seq file, seed.fa, must exist."
    print >> sys.stderr, "  This script, the hmm file, and the seq file"
    print >> sys.stderr, "  must all exist in the same directory."
    sys.exit(0)

  hmm_filename = sys.argv[1]
  hmm_basename = os.path.basename(hmm_filename)
  seed = os.path.splitext(hmm_basename)[0]
  runname = "pdb_against_%s_hmm" % (seed)

  # Compare HMM against all PDB sequences
  out_filename = "%s.out" % (runname)
  out_file = open(out_filename, "w")
  print >>  out_file, "Comparing HMM against PDB"
  cmd = "hmmscore %s -dbsize 100000 -db %s -sw 2 -i %s 2>> %s" %  \
      (runname, pdb_db_filename, hmm_filename, out_filename)
  print >> out_file, "Running the following command:"
  print >> out_file, cmd
  print >> out_file, "hmmscore output:"
  out_file.close()
  os.system(cmd)
  out_file = open(out_filename, "a")

  # Extract significant hits from output
  print >>  out_file, "Extracting significant hits from hmmscore output."
  pdb_hits_filename = "%s.dist" % (runname)
  if not os.path.exists(pdb_hits_filename):
    print >> out_file, "hmmscore produced no output file"
    sys.exit(0)
  hmmscore_output_filesize = os.path.getsize(pdb_hits_filename)
  if (hmmscore_output_filesize) == 0:
    print >> out_file, "%s  is empty" % (pdb_hits_filename)
  sig_hit_ids =  get_significant_hits_from_hmmscore_output(pdb_hits_filename,
     eval_threshold)
  fastacmd_in_file = open(fastacmd_in_filename, "w")
  num_sig_hits = len(sig_hit_ids)
  print >> out_file, "%d hits found in %d with eval exponent < %d" % \
    (num_sig_hits, pdb_db_filename, eval_threshold)
  if num_sig_hits == 0:
    print >> sys.stderr, "No hits with eval exponent < %d" % (eval_threshold)
    sys.exit(0)
  for hit in sig_hit_ids:
    print >> fastacmd_in_file, hit
  fastacmd_in_file.close()

  # Run fastacmd to get pdb sequences for top hits
  print >> out_file, "Retrieving PDB sequences for top hits"
  print >> out_file, "Running the following command:"
  cmd = "fastacmd -i %s -d %s >| %s" % (fastacmd_in_filename,
     pdb_db_filename, fastacmd_out_filename)
  print >> out_file, cmd
  os.system(cmd)
  #os.system("rm %s" % (fastacmd_in_filename))

  # Run align2model, including seed in alignment
  print >> out_file, "Aligning top hits to input HMM"
  print >> out_file, "Running the following command:"
  cmd = "align2model %s -sw 2 -i %s -db %s -db %s.fa -adpstyle 5 >& %s.align2model.log" % \
      (runname, hmm_filename, fastacmd_out_filename, seed, runname)
  print >> out_file, cmd
  os.system(cmd)
  out_file.close()

if __name__ == "__main__":
  main()
