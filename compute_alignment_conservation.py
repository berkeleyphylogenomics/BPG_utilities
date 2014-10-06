#!/usr/bin/env python

import os, re, string
from Bio import AlignIO, SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.pairwise2 import dictionary_match
from optparse import OptionParser
from pfacts003.phylofacts.models import Family, TreeNodeAlignmentConservation

blosum62_of_residues = dictionary_match(SubsMat.SeqMat(MatrixInfo.blosum62))

def get_alignment_seqs_and_aligned_column_indices(alignment):
  alignment_length = 0
  aligned_column_indices = set()
  alignment_seqs = {}
  first_pass = True
  i = 0

  for row in alignment:
    seq = row.seq.tostring()
    if first_pass:
      alignment_length = len(row.seq)
      for j in range(len(seq)):
        if seq[j] == '-' or seq[j].isupper():
          aligned_column_indices.add(j)
      first_pass = False
    alignment_seqs[i] = seq
    i += 1
  return (alignment_seqs, aligned_column_indices)

def get_conservation_info(alignment_seqs, aligned_column_indices):
  column_conserved_residue = {}
  column_score = {}
  for j in aligned_column_indices:
    freq_of_residue = {}
    highest_frequency = 0
    most_frequent_residue = ''
    for i in alignment_seqs.keys():
      residue = alignment_seqs[i][j]
      if residue == '-':
        continue
      if residue not in freq_of_residue:
        freq_of_residue[residue] = 0
      freq_of_residue[residue] += 1
      if freq_of_residue[residue] > highest_frequency:
        highest_frequency = freq_of_residue[residue]
        most_frequent_residue = residue
    column_conserved_residue[j] = most_frequent_residue
    num_pairs = 0
    sum_of_scores = 0.0
    for i0 in range(len(alignment_seqs)):
      residue0 = alignment_seqs[i0][j]
      if residue0 != '-':
        for i1 in range(i0):
          residue1 = alignment_seqs[i1][j]
          if residue1 != '-':
            score = blosum62_of_residues(alignment_seqs[i0][j],
                                          alignment_seqs[i1][j])
            sum_of_scores += score
            num_pairs += 1
    if num_pairs > 0:
      column_score[j] = sum_of_scores / num_pairs
    else:
      print "Gappy column %d" % j
  return (column_conserved_residue, column_score)

def update_tree_node_alignment_conservation(tree_node, aligned_column_indices,
                                          column_conserved_residue,
                                          column_score):
  for column_index in aligned_column_indices:
    objects = TreeNodeAlignmentConservation.objects.filter(
                            tree_node = tree_node, column_index = column_index)
    if objects:
      object = objects[0]
      object.conserved_residue = column_conserved_residue[column_index]
      object.blosum62_conservation_score = column_score[column_index]
      object.save()
    else:
      TreeNodeAlignmentConservation.objects.create(tree_node = tree_node,
                    column_index = column_index,
                    conserved_residue = column_conserved_residue[column_index],
                    blosum62_conservation_score = column_score[column_index])

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if len(args[0]) < 10 or args[0][0:3] != 'bpg':
    opt_parser.error('Argument must be a bpg accession like bpg0123456')
  bpg_accession = args[0]
  try:
    family_id = int(bpg_accession[3:])
  except ValueError:
    opt_parser.error('Argument must be a bpg accession like bpg0123456')
  family_dir = '/clusterfs/ohana/bpg/pfacts/%s/%s/%s' % (bpg_accession[0:4],
                                                          bpg_accession[0:7],
                                                          bpg_accession) 
  if not os.path.exists(family_dir):
    opt_parser.error('Family %s not found on the filesystem.' % bpg_accession)
  family = Family.objects.get(id = family_id)
  if family.status == 'bad':
    opt_parser.error('Family %s is marked as bad in the database.'  \
                      % bpg_accession)

  alignment_file = os.path.join(family_dir, '%s.a2m' % bpg_accession)
  f = open(alignment_file)
  for alignment in AlignIO.parse(f, "fasta"):
    break

  alignment_seqs, aligned_column_indices \
      = get_alignment_seqs_and_aligned_column_indices(alignment)
  column_conserved_residue, column_score \
      = get_conservation_info(alignment_seqs, aligned_column_indices)

  outfname = os.path.join(family_dir,
                        bpg_accession + '.alignmentconservation.csv')
  outf = open(outfname, 'w')
  outf.write('ColumnIndex,ConservedResidue,Blosum62ConservationScore\n')
  for j in aligned_column_indices:
    outf.write('%d,%s,%f\n' % (j, column_conserved_residue[j], column_score[j]))
  outf.close()

  root = family.canonical_root_node()
  update_tree_node_alignment_conservation(root, aligned_column_indices,
                                          column_conserved_residue,
                                          column_score)

if __name__ == '__main__':
  main()
