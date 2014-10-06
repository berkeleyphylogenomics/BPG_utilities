#!/usr/bin/env python

import glob
import os
import string
import sys

from Bio import SeqIO
from Bio.SeqUtils import CheckSum

from bpg.common import BPGPWID
from pfacts003.phylofacts.models import Sequence, SequenceHeader, UniProt, \
TreeNode, TreeNodeAlignment, Family

def main():
  basepath = '/clusterfs/ohana/bpg/Hpylori26695/GHGs'
  os.chdir(basepath)
  f = open('Helicobacter_pylori_26695')
  input_records = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
  f.close()
  set_cover_dict = {}
  seed_dict = {}
  f = open("clusters")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line.split()[0] == 'Cluster':
      starting_cluster = True
    else:
      if starting_cluster:
        seed_id = line.split()[0]
        set_cover_dict[seed_id] = seed_id
        seed_dict[seed_id] = set([seed_id])
        starting_cluster = False
      else:
        seq_id = line.split()[0]
        set_cover_dict[seq_id] = seed_id
        seed_dict[seed_id].add(seq_id)
  uppercase_translation = string.maketrans(string.lowercase, string.uppercase)
  dotdash='.-'
  seed_paths = glob.glob('seeds*/*')
  for seed_path in seed_paths:
    seed_dir, seed = os.path.split(seed_path)
    sought_seguids = [(id, CheckSum.seguid(input_records[id].seq))
                      for id in
                      seed_dict[('lcl|%s' % seed)]]
    alignment_paths = glob.glob('%s/bpg*.a2m' % seed_path)
    if alignment_paths:
      family_id = int(os.path.splitext(os.path.split(alignment_paths[0])[1])[0][3:])
      family = Family.objects.get(id = family_id)
      root = family.canonical_root_node()
      tree_node_objects = TreeNode.objects.filter(tree = family.canonical_tree,
                            sequence_header__isnull = False).order_by('left_id')
      sequence_headers = [node.sequence_header for node in tree_node_objects]
      tree_node_alignment_objs = \
        TreeNodeAlignment.objects.filter(tree_node = root)
      alignment_of_sequence_header = {}
      for obj in tree_node_alignment_objs:
        alignment_of_sequence_header[obj.sequence_header] = obj
      sequence_headers_of_seguid = {}
      for obj in tree_node_alignment_objs:
        seguid = obj.sequence_header.sequence.seguid
        if seguid not in sequence_headers_of_seguid:
          sequence_headers_of_seguid[seguid] = set()
        sequence_headers_of_seguid[seguid].add(obj.sequence_header)
      for id, seguid in sought_seguids:
        if seguid in sequence_headers_of_seguid:
          uniprot_identifiers = list(set(
                    [sequence_header.uniprot.uniprot_identifier for
                    sequence_header in sequence_headers_of_seguid[seguid]
                    if sequence_header.uniprot]))
          if uniprot_identifiers:
            print '%s,"%s"' % (id, ','.join(uniprot_identifiers))
          else:
            the_sequence_header = list(sequence_headers_of_seguid[seguid])[0]
            seq0 = alignment_of_sequence_header[
                                    the_sequence_header].aligned_sequence.chars
            max_percent_id = 0.0
            closest_uniprot = None
            for obj in tree_node_alignment_objs:
              if obj.sequence_header.uniprot:
                seq1 = obj.aligned_sequence.chars
                pwid = BPGPWID.pairwise_identity_KS_1(seq0, seq1)
                if pwid > max_percent_id:
                  max_percent_id = pwid
                  closest_uniprot = obj.sequence_header.uniprot
            print '%s,"%s(%0.3f)"' % (id, closest_uniprot.uniprot_identifier,
                                      max_percent_id)
        else:
          print "No exact match for %s" % id
    else:
      alignment_path = '%s/final.a2m' % seed_path
      if os.path.exists(alignment_path):
        uniprot_identifiers_of_seguid = {}
        f = open(alignment_path)
        flowerpower_alignment = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        f.close()
        for id in flowerpower_alignment:
          if id[0:3] == 'tr|' or id[0:3] == 'sp|':
            unaligned_seq = flowerpower_alignment[
                id].seq.tostring().translate(uppercase_translation, dotdash)
            seguid = CheckSum.seguid(unaligned_seq)
            if seguid not in uniprot_identifiers_of_seguid:
              uniprot_identifiers_of_seguid[seguid] = set()
            uniprot_identifiers_of_seguid[seguid].add(id.split('|')[2])
        for id, seguid in sought_seguids:
          if seguid in uniprot_identifiers_of_seguid:
            print '%s,"%s"' % (id,
                ','.join(list(uniprot_identifiers_of_seguid[seguid])))
          else:
            print "No exact match for %s" % id
        
      else:
        print "Seed %s has not been FlowerPowered" % seed

if __name__ == '__main__':
  main()
