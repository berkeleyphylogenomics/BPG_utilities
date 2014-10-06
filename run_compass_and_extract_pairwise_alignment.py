#!/usr/bin/env python

import commands, sys
from array import array
from optparse import OptionParser
from Bio import Seq, SeqIO, AlignIO

def main():
  usage = "%prog [options] sequence_id_1 msa_1 sequence_id_2 msa_2"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-v", "--verbose", dest="verbose",
              action="store_true", default=True,
              help="Whether to print verbose output")
  opt_parser.add_option("-q", "--quiet", dest="verbose",
              action="store_false", default=True,
              help="Whether to suppress verbose output")
  (options, args) = opt_parser.parse_args()
  if len(args) != 4:
    opt_parser.error('Incorrect number of arguments')
  seed1 = args[0]
  msa1 = args[1]
  seed2 = args[2]
  msa2 = args[3]
  cmd = "compass -c 0 -i %s -j %s" % (msa1, msa2)
  status, output = commands.getstatusoutput(cmd)
  f = open("%s_%s_compass.out" % (seed1, seed2), "w")
  f.write(output)
  f.close()
  lines = output.split('\n')
  seeking_alignment = True
  used_seq_id1 = ''
  used_seq_id2 = ''
  which_seq = 1
  aligned_seq1 = array('c')
  aligned_seq2 = array('c')
  aligned_seq1_start_index = -1
  aligned_seq2_start_index = -1
  aligned_seq1_prev_index = -1
  aligned_seq2_prev_index = -1
  for line in lines:
    if seeking_alignment:
      if line[0:5] == 'Smith':
        seeking_alignment = False
    else:
      if line[0:10] == 'Parameters':
        break
      else:
        line = line.rstrip()
        if len(line) > 0:
          if len(line[0:23].rstrip()) > 0:
            if which_seq == 1:
              if used_seq_id1 == '':
                used_seq_id1 = line[0:23].rstrip()
                aligned_seq1_start_index = int(line[23:30].rstrip()) - 1
              seq, cur_end_index = line[30:].split()
              aligned_seq1.fromstring(seq)
              cur_start_index = int(line[23:30].rstrip())
              if aligned_seq1_prev_index > 0 and \
                   cur_start_index != aligned_seq1_prev_index + 1:
                print "Alert! Alignment break in seq1 at %d,%d" \
                    % (aligned_seq1_prev_index, cur_start_index)
              aligned_seq1_prev_index = int(cur_end_index)
              which_seq = 2
            else:
              if used_seq_id2 == '':
                used_seq_id2 = line[0:23].rstrip()
                aligned_seq2_start_index = int(line[23:30].rstrip()) - 1
              seq, cur_end_index = line[30:].split()
              aligned_seq2.fromstring(seq)
              cur_start_index = int(line[23:30].rstrip())
              if aligned_seq2_prev_index > 0 and \
                   cur_start_index != aligned_seq2_prev_index + 1:
                print "Alert! Alignment break in seq2 at %d,%d" \
                    % (aligned_seq2_prev_index, cur_start_index)
              aligned_seq2_prev_index = int(cur_end_index)
              which_seq = 1
  if len(aligned_seq1) != len(aligned_seq2):
    print "Error: lengths of COMPASS aligned sequences differ"
    sys.exit(1)
  f = open(msa1)
  alignments1 = list(AlignIO.parse(f, "fasta"))
  alignment1 = SeqIO.to_dict(alignments1[0])
  f.close()
  f = open(msa2)
  alignments2 = list(AlignIO.parse(f, "fasta"))
  alignment2 = SeqIO.to_dict(alignments2[0])
  f.close()
  aligned_seed_seq1 = alignment1[seed1].seq.tostring()
  input_aln_of_used_seq1 = alignment1[used_seq_id1].seq.tostring()
  aligned_seed_seq2 = alignment2[seed2].seq.tostring()
  input_aln_of_used_seq2 = alignment2[used_seq_id2].seq.tostring()
  i = 0
  j = 0
  k = 0
  l = 0
  new_aligned_seed_seq1 = array('c')
  new_aligned_seed_seq2 = array('c')
  while j < aligned_seq1_start_index:
    if aligned_seed_seq1[i].isalpha():
      new_aligned_seed_seq1.fromstring(aligned_seed_seq1[i].upper())
      new_aligned_seed_seq2.fromstring('-')
    if input_aln_of_used_seq1[i].isalpha():
      j += 1
    i += 1
  while not input_aln_of_used_seq1[i].isalpha():
    if aligned_seed_seq1[i].isalpha():
      new_aligned_seed_seq1.fromstring(aligned_seed_seq1[i].upper())
      new_aligned_seed_seq2.fromstring('-')
    i += 1
  while k < aligned_seq2_start_index:
    if aligned_seed_seq2[l].isalpha():
      new_aligned_seed_seq1.fromstring('-')
      new_aligned_seed_seq2.fromstring(aligned_seed_seq2[l].upper())
    if input_aln_of_used_seq2[l].isalpha():
      k += 1
    l += 1
  while not input_aln_of_used_seq2[l].isalpha():
    if aligned_seed_seq2[l].isalpha():
      new_aligned_seed_seq1.fromstring('-')
      new_aligned_seed_seq2.fromstring(aligned_seed_seq2[l].upper())
    l += 1
  for n in xrange(len(aligned_seq1)):
    if aligned_seq1[n].isalpha() or aligned_seq2[n].isalpha():
      if aligned_seq1[n].isalpha():
        if aligned_seq1[n].upper() != input_aln_of_used_seq1[i].upper():
          print 'Error! Mismatch at COMPASS col %d, MSA1 col %d' % (n, i)
          sys.exit(1)
        new_aligned_seed_seq1.fromstring(aligned_seed_seq1[i].upper())
        i += 1
      else:
        new_aligned_seed_seq1.fromstring('-')
      if aligned_seq2[n].isalpha():
        if aligned_seq2[n].upper() != input_aln_of_used_seq2[l].upper():
          print 'Error! Mismatch at COMPASS col %d, MSA2 col %d' % (n, l)
          sys.exit(1)
        new_aligned_seed_seq2.fromstring(aligned_seed_seq2[l].upper())
        l += 1
      else:
        new_aligned_seed_seq2.fromstring('-')
      while i < len(aligned_seed_seq1) and \
          not input_aln_of_used_seq1[i].isalpha():
        if aligned_seed_seq1[i].isalpha():
          new_aligned_seed_seq1.fromstring(aligned_seed_seq1[i].upper())
          new_aligned_seed_seq2.fromstring('-')
        i += 1
      while l < len(aligned_seed_seq2) and \
          not input_aln_of_used_seq2[l].isalpha():
        if aligned_seed_seq2[l].isalpha():
          new_aligned_seed_seq1.fromstring('-')
          new_aligned_seed_seq2.fromstring(aligned_seed_seq2[l].upper())
        l += 1
  while i < len(aligned_seed_seq1):
    if aligned_seed_seq1[i].isalpha():
      new_aligned_seed_seq1.fromstring(aligned_seed_seq1[i].upper())
      new_aligned_seed_seq2.fromstring('-')
    i += 1
  while l < len(aligned_seed_seq2):
    if aligned_seed_seq2[l].isalpha():
      new_aligned_seed_seq1.fromstring('-')
      new_aligned_seed_seq2.fromstring(aligned_seed_seq2[l].upper())
    l += 1
  print ">%s" % seed1
  print new_aligned_seed_seq1.tostring()
  print ">%s" % seed2
  print new_aligned_seed_seq2.tostring()


if __name__ == '__main__':
  main()
