#!/usr/bin/env python

import os, subprocess, cPickle


def main():
  f = open('/clusterfs/vasudha/bpg/OrthologsForQuest/threshold_of_taxon_pair.pkl')
  threshold_of_taxon_pair = cPickle.load(f)
  f.close()
  taxa = threshold_of_taxon_pair.keys()
  taxon_id1 = taxa[int(os.environ['PBS_ARRAYID'])]
  os.chdir('/clusterfs/vasudha/bpg/OrthologsForQuest/orthoXML')
  for taxon_id2 in threshold_of_taxon_pair[taxon_id1]:
    p = subprocess.Popen(['write_pairwise_species_file.py',
                          str(taxon_id1), str(taxon_id2)])
    status = os.waitpid(p.pid, 0)[1]

if __name__ == '__main__':
  main()
