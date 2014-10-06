#!/usr/bin/env python

import os, sys, cPickle
from Bio import SeqIO

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <flock_dir> <fasta_of_genes_to_cover>" % sys.argv[0]
    sys.exit(0)

  flock_dir = sys.argv[1]
  f = open(os.path.join(flock_dir, "set_cover_dict.pkl"))
  set_cover_dict = cPickle.load(f)
  f.close()
  f = open(os.path.join(flock_dir, "cluster_dict.pkl"))
  cluster_dict = cPickle.load(f)
  f.close()

  f = open(sys.argv[2])
  records_to_cover = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
  f.close()

  accessions_to_cover = set()
  target_accessions = set(['lcl|%s' % id.replace('|','_') for 
                          id in records_to_cover.keys()])
  target_genes_covered_by_seed = {}

  max_num_target_genes_covered_by_seed = 0
  num_genes_covered_by_max_covering_seed = 0
  max_covering_seed = ""
  max_covering_seed_is_target = False
  flock_seed_is_target = {}
  print "Missing accessions"
  for accession in target_accessions:
    if accession in set_cover_dict:
      accessions_to_cover.add(accession)
      for flock_seed in set_cover_dict[accession]:
        if flock_seed not in target_genes_covered_by_seed:
          target_genes_covered_by_seed[flock_seed] = set()
        target_genes_covered_by_seed[flock_seed].add(accession)
        num_covered = len(target_genes_covered_by_seed[flock_seed])
        flock_seed_is_target[flock_seed] = (flock_seed in target_accessions)
        if num_covered > max_num_target_genes_covered_by_seed or \
            num_covered == max_num_target_genes_covered_by_seed and \
            (not max_covering_seed_is_target and \
                flock_seed_is_target[flock_seed] \
            or len(cluster_dict[flock_seed]) \
                  > num_genes_covered_by_max_covering_seed):
          max_num_target_genes_covered_by_seed = num_covered
          max_covering_seed = flock_seed
          max_covering_seed_is_target = flock_seed_is_target[flock_seed]
          num_genes_covered_by_max_covering_seed = len(cluster_dict[flock_seed])
    else:
      print "%s,???" % accession

  covered_accessions = set()
  seeds_used = set()
  print "Covering seeds"
  while len(accessions_to_cover) > 0 and \
        max_num_target_genes_covered_by_seed > 0:
    print "%s,%d,%d" % (max_covering_seed, max_num_target_genes_covered_by_seed,
                        num_genes_covered_by_max_covering_seed)
    covered_accessions |= target_genes_covered_by_seed[max_covering_seed]
    accessions_to_cover -= target_genes_covered_by_seed[max_covering_seed]
    seeds_used.add(max_covering_seed)
    max_num_target_genes_covered_by_seed = 0
    max_covering_seed = ""
    max_covering_seed_is_target = False
    num_genes_covered_by_max_covering_seed = 0
    for seed in target_genes_covered_by_seed.keys():
      if seed not in seeds_used:
        newly_covered = target_genes_covered_by_seed[seed] & accessions_to_cover
        if len(newly_covered) > max_num_target_genes_covered_by_seed or \
            len(newly_covered) == max_num_target_genes_covered_by_seed and \
            (len(cluster_dict[seed]) > num_genes_covered_by_max_covering_seed \
              or \
             len(cluster_dict[seed]) == num_genes_covered_by_max_covering_seed \
             and not max_covering_seed_is_target \
             and flock_seed_is_target[seed]):
          max_num_target_genes_covered_by_seed = len(newly_covered)
          max_covering_seed = seed
          max_covering_seed_is_target = flock_seed_is_target[seed]
          num_genes_covered_by_max_covering_seed = len(cluster_dict[seed])

  if len(accessions_to_cover) > 0:
    print "%d uncovered accessions" % len(accessions_to_cover)
    for accession in accessions_to_cover:
      print accession

if __name__ == '__main__':
  main()
