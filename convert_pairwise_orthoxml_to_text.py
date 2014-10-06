#!/usr/bin/env python2.6

import sys
import xml.etree.cElementTree as ET
from pfacts003.phylofacts.models import UniProt

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <orthoXML_file> <output_file>" % sys.argv[0]
    sys.exit(0)
  outf = open(sys.argv[2], 'w')
  xmltree = ET.parse(sys.argv[1])
  uniprot_accession_of_id = {}
  ns = '{http://orthoXML.org/2011/}'
  num_gene_nodes = 0
  for geneNode in xmltree.getiterator('%sgene' % ns):
    uniprot_id = geneNode.get('id')
    uniprot_accession = geneNode.get('protId')
    uniprot_accession_of_id[uniprot_id] = uniprot_accession

  outf.write('GeneAccession,GeneDescription,OrthologAccession,Score,Support\n')
  for orthologGroupNode in xmltree.getiterator('%sorthologGroup' % ns):
    first = True
    for geneRef in orthologGroupNode:
      uniprot_id = geneRef.get('id')
      if first:
        uniprot = UniProt.objects.get(id = uniprot_id)
        outf.write('%s,%s,' % (uniprot.accession, uniprot.de))
        first = False
      else:
        uniprot_accession = uniprot_accession_of_id[uniprot_id]
        outf.write('%s,' % uniprot_accession)
        for score in geneRef.getiterator('%sscore' % ns):
          outf.write('%s,' % score.get('value'))
        outf.write('"')
        for supporting_phog in geneRef.getiterator('%ssupporting_phog' % ns):
          outf.write('%s:%s:%s; ' % (supporting_phog.get('family'),
                                    supporting_phog.get('phog'),
                                    supporting_phog.get('minimum_threshold')))
        outf.write('"')
    outf.write('\n')
  outf.close()
      

if __name__ == '__main__':
  main()
