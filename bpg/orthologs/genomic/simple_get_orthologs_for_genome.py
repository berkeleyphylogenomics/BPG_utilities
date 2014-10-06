#!/usr/bin/env python

import os, sys, glob
from pfacts003.phylofacts.models import *
from pfacts003.phog.orthologs import getOrthologQuerySet

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <taxonomy_id_of_genome>" % sys.argv[0]
    sys.exit(0)

  taxonomy_id = sys.argv[1]

  dir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5'
  os.chdir(dir)
  genome_files = glob.glob('%s_*.fasta' % taxonomy_id)
  if len(genome_files) == 0:
    print "No genome file found for %s" % taxonomy_id
    sys.exit(1)
  genome_file = genome_files[0]

  f = open(genome_files[0])
  accessions = [line.strip().split()[0].split(':')[1] for line in f.readlines()
                if len(line) > 0 and line[0] == '>']
  f.close()

  taxon_ids = [ 83333, 35554, 289376, 63363, 1773, 1148, 33072, 324602,
  1423, 2336, 515635, 1299, 478009, 2214, 2190, 311400, 2287, 374847 ]
  interesting_taxon_ids = set(taxon_ids)
  sys.stdout.write('UniProtAccession,Description,')
  sys.stdout.write(','.join(['OrthologsInTaxon%d' % taxon_id for taxon_id in
                                                                taxon_ids]))
  sys.stdout.write('\n')
  for accession in accessions:
    uniprot_dat_indices = UniProtDatIndex.objects.filter(uniprot_accession =
                                                            accession)
    if uniprot_dat_indices:
      orthologs_in_taxon = {}
      for taxon_id in interesting_taxon_ids:
        orthologs_in_taxon[taxon_id] = set()
      tree_nodes = TreeNode.objects.filter(sequence_header__uniprot
                                  = uniprot_dat_indices[0].uniprot).filter(
                                  tree__method = 'ml').exclude(
                                  tree__family__status__exact = 'bad')
      for tree_node in tree_nodes:
        containing_phog = tree_node.get_containing_phog(threshold = 0.125)
        if containing_phog:
          interesting_orthologs = containing_phog.get_included_leaves().filter(
                sequence_header__uniprot__taxon__id__in = interesting_taxon_ids)
          for ortholog in interesting_orthologs:
            orthologs_in_taxon[ortholog.sequence_header.uniprot.taxon.id].add(
                                ortholog.sequence_header.uniprot)
        
      sys.stdout.write('%s,"%s",' % (accession, 
                                      uniprot_dat_indices[0].uniprot.de))
      sys.stdout.write(','.join(['"%s"' % ';'.join(
          [ortholog.accession for ortholog in orthologs_in_taxon[taxon_id]]) 
            for taxon_id in taxon_ids]))
      sys.stdout.write('\n')

if __name__ == '__main__':
  main()
