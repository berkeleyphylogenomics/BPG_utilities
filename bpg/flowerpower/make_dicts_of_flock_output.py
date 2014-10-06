#!/usr/bin/env python

import cPickle
import os, glob
from Bio import SeqIO

def main():
  root_dir = os.getcwd()
  seeds_dirs = glob.glob('seeds*')
  seed_dir_dict = {}
  set_cover_dict = {}
  cluster_dict = {}
  for seeds_dir in seeds_dirs:
    os.chdir(seeds_dir)
    seeds = os.listdir('.')
    for seed in seeds:
      seed_dir_dict[seed] = '%s/%s/%s' % (root_dir, seeds_dir, seed)
      cluster_dict[seed] = set()
      cluster_file = '%s/%s_cluster.fa' % (seed, seed)
      if os.path.exists(cluster_file):
        f = open(cluster_file, "r")
        seq_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        f.close()
        for id in seq_dict.keys():
          cluster_dict[seed].add(id)
          if id not in set_cover_dict:
            set_cover_dict[id] = set()
          set_cover_dict[id].add(seed)
      else:
        cluster_dict[seed].add(seed)
        if seed not in set_cover_dict:
          set_cover_dict[seed] = set()
        set_cover_dict[seed].add(seed)
    os.chdir('..')
  print "Done gathering dictionaries, now pickling them..."
  f = open("seed_dir_dict.pkl", "w")
  cPickle.dump(seed_dir_dict, f)
  f.close()
  f = open("cluster_dict.pkl", "w")
  cPickle.dump(cluster_dict, f)
  f.close()
  f = open("set_cover_dict.pkl", "w")
  cPickle.dump(set_cover_dict, f)
  f.close()

if __name__ == '__main__':
  main()
