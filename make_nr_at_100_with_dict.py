#!/usr/bin/python

import os, sys
from Bio import SeqIO
from Bio.SeqUtils import CheckSum
import cPickle

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <runname> <fasta_file>" % sys.argv[0]
    sys.exit(0)
  runname = sys.argv[1]
  input_fasta = sys.argv[2]
  input_records = SeqIO.parse(open(input_fasta, "rU"), "fasta")
  output_handle = open("%s_nr100.fa" % runname, "w")
  nr_seqs={}
  dups_of_seq = {}
  for input_record in input_records:
    header = input_record.id
    if header.find('|') >= 0:
      header = header.split('|')[1]
    header = header.replace("%", "%25")
    header = header.replace(":", "%3a")
    header = header.replace(";", "%3b")
    header = header.replace(",", "%2c")
    header = header.replace("(", "%28")
    header = header.replace(")", "%29")
    the_seguid = CheckSum.seguid(input_record.seq)
    seqstr = input_record.seq.tostring()
    if the_seguid in nr_seqs and nr_seqs[the_seguid] != None:
      if seqstr in nr_seqs[the_seguid] \
          and nr_seqs[the_seguid][seqstr] != None:
        nr_seqs[the_seguid][seqstr].add(header)
        dups_of_seq[header] = nr_seqs[the_seguid][seqstr]
      else:
        nr_seqs[the_seguid][seqstr] = set([header])
        SeqIO.write([input_record], output_handle, "fasta")
        dups_of_seq[header] = nr_seqs[the_seguid][seqstr]
    else:
      nr_seqs[the_seguid] = { seqstr: set([header]) }
      SeqIO.write([input_record], output_handle, "fasta")
      dups_of_seq[header] = nr_seqs[the_seguid][seqstr]

  pklfp = open("%s_dict.pkl" % runname, "w")
  cPickle.dump((dups_of_seq, nr_seqs), pklfp)
  pklfp.close()
        
if __name__ == '__main__':
  main()
