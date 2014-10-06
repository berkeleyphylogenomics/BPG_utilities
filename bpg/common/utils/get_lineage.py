#!/usr/bin/python

import glob, sys, os.path
from pfacts003.phylofacts.models import *

'''
This script accepts a taxon id, finds the file with the genome for that
taxon in the /clusterfs/ohana/external/genomes/QuestForOrthologs/Release5/
directory, figures out if each accession is found in any family, and
prints out  a report on the coverage in PhyloFacts for the genome.
'''

if  len(sys.argv) < 2:
    print "usage: specify taxonomy id"
    sys.exit(1)

taxon_id = sys.argv[1]
       
taxon = UniProtTaxonomy.objects.get(id = taxon_id)

phylum=UniProtTaxonomy.objects.filter (left_id__lte = taxon.left_id,right_id__gte = taxon.right_id,rank ='Phylum')

kingdom=UniProtTaxonomy.objects.filter (left_id__lte = taxon.left_id,right_id__gte = taxon.right_id,rank ='Kingdom')

#lineage = taxon.lineage()
#print "Taxonomy: %s" % "/".join([str(item) for item in lineage])

print "Species: %s, Phylum: %s" % (taxon.scientific_name, Kingdom)

