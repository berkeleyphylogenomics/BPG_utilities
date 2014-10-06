#!/usr/bin/python

import sys
from pdb2scop import *

"""
Input: a concatenation of BLAST outputs for PhyloFacts ghmm
  consensus sequences vs. Astral PDB40.
Output: a graphical representation of how well covered
  each query sequence is by its BLAST hits
"""

pdb40_filename = "/home/ruchira/CASP8/astral40.fasta"
fastacmd_out_filename = "fastacmd.out"
hmmscore_run_name = "tmp"
hmmscore_out_filename = "%s.dist" % (hmmscore_run_name)

def get_blast_start_end(hit_list, subject):
  blast_query_start = None
  blast_query_end = None
  #print  "subject =", subject
  for hit in hit_list:
    id = hit[0]
    #print "id =", id
    if id  == subject:
      blast_query_start = hit[1]
      blast_query_end = hit[2]
      break;
  return blast_query_start, blast_query_end 

def parse_eval(eval_string):
  if eval_string.find("e") > -1:
    eval = eval_string.split("e")
    coefficient = float(eval[0])
    exponent = int(eval[1])
  else:
    coefficient = float(eval_string)
    exponent = 0
  return coefficient, exponent

def main():

  # Open file and read one line at a time.
  infile = open(sys.argv[1], "r")
  line = infile.readline()
  last_query = ""
  new_query = True
  
  while line != "":
    # pass over comment lines except for one displaying query length
    while line != "" and line[0] == "#":
      #print line
      # query_seqlen is in a line that begins with "#len"
      # read it when found
      if line[1:4] == "len":
        query_seqlen = int(line[4:])
      line = infile.readline()
    
    # set up stuff for this query
    pdb_id_filename = "pdb_id.list"
    pdb_id_uniq_filename = "pdb_id_uniq.list"
    pdb_id_file = open(pdb_id_filename, "w")
    hit_list = []
    
    # store the blast results for this query
    while line != "" and line[0] != "#":
      #print line
      # split line into tokens by tabs
      tokens = line.split("\t")
      this_query = tokens[0]
      subject = tokens[1]
      blast_query_start = int(tokens[6])
      blast_query_end = int(tokens[7])
      blast_eval = tokens[10]
      hit_list.append((subject, blast_query_start, blast_query_end, blast_eval))
      print >> pdb_id_file, subject
      line = infile.readline()

    print >> sys.stdout, "#====== %s len %d ======" % (this_query, query_seqlen)
    
    pdb_id_file.close()
    cmd = "sort %s | uniq > %s" % (pdb_id_filename, pdb_id_uniq_filename)
    os.system(cmd)

    # process these blast results
    # get the sequences for all the subjects
    cmd = "fastacmd -i %s -d %s > %s" % (pdb_id_uniq_filename, pdb40_filename,
         fastacmd_out_filename)
    #print cmd
    os.system(cmd)
    
    # score against hmm for query BPG family
    query_data_dir = "/home/bpg/pfacts/%s" % (this_query)
    if os.path.exists("%s/user" % (query_data_dir)):
      query_hmm_subdir = "user"
    elif os.path.exists("%s/T2K" % (query_data_dir)):
      query_hmm_subdir = "T2K"
    elif os.path.exists("%s/FP" % (query_data_dir)):
      query_hmm_subdir = "FP"
    else:
      # skip this query; we can't find an hmm
      print >> sys.stderr, "no user, T2K, or FP subdir for %s" % (this_query)
      continue
    
    # *** FIRST, CALIBRATE HMM -- how to do?
    cmd = "hmmscore %s  -i %s/%s/%s.mod -db %s -sw 2 -dbsize 100000 &> hmmscore.out" % \
        (hmmscore_run_name, query_data_dir, query_hmm_subdir, this_query, fastacmd_out_filename)
    #print cmd
    os.system(cmd)
    
    # read hmmscore output file and output info on good hits
    hmmscore_out_file = open(hmmscore_out_filename, "r")
    hmm_output_line = hmmscore_out_file.readline()
    # skip comment lines
    while hmm_output_line != "" and hmm_output_line[0] == "%":
      hmm_output_line = hmmscore_out_file.readline()
    while hmm_output_line != "":
      elements = hmm_output_line.split()
      [subject, length, simple, reverse, hmm_eval, domain] = \
          hmm_output_line.split()[:6]
      subject = subject.split("|")[1]  #take off prefix
      domain = domain[1:]   #take off leading colon
      hmmscore_eval_coefficient, hmmscore_eval_exponent = parse_eval(hmm_eval)
      if hmmscore_eval_exponent < -2:
        blast_query_start, blast_query_end = get_blast_start_end(hit_list, subject)
        if blast_query_start == None or blast_query_end == None:
          # skip this subject
          print >> sys.stderr, "Missing blast start and/or end position."
        else:
          # print subject, domain, eval, blast_query_start, blast_query_end
          sys.stdout.write("%-9s %-11s %de%4d %5d %5d " % (subject, domain,
                  hmmscore_eval_coefficient, hmmscore_eval_exponent,
                  blast_query_start, blast_query_end))
          # print line of dashes showing coverage of this hit
          for i in range (1,blast_query_start):
            if i%10==0: sys.stdout.write(" ")
          for i in range (blast_query_start,blast_query_end):
            if i%10==0: sys.stdout.write("-")
          sys.stdout.write("\n")

      hmm_output_line = hmmscore_out_file.readline()
  
    last_query = this_query

  
  
if __name__ == "__main__":
  main()
