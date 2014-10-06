#!/usr/bin/env python

import sys, cPickle
from pfacts003.phylofacts.models import UniProt, UniProtTaxonomy, TreeNode

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <idmap_file>" % sys.argv[0]
    sys.exit(0)

  f = open(sys.argv[1])
  idmap = cPickle.load(f)
  f.close()
  headers = idmap.values()
  headers_of_taxon = {}
  for header in headers:
    id = header.split()[0]
    id_components = id.split('|')
    if len(id_components) == 3:
      uniprot_accession = id_components[1]
      uniprot = UniProt.objects.get(accession = uniprot_accession)
      if uniprot.taxon not in headers_of_taxon:
        headers_of_taxon[uniprot.taxon] = set()
      headers_of_taxon[uniprot.taxon].add(header)

  ordered_taxa = [(taxon.left_id, taxon) 
                  for taxon in headers_of_taxon.keys()]
  ordered_taxa.sort()
  padding = ''.join([' ' for x in range(40)])
  for left_id, taxon in ordered_taxa:
    print ':'.join([t.__unicode__() for t in taxon.lineage()])
    print
    for header in headers_of_taxon[taxon]:
      print "%s%s" % (padding, header)
    print

if __name__ == '__main__':
  main()
