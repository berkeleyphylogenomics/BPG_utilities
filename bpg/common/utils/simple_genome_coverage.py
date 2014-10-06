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
files = glob.glob("/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5/%s_*.fasta" 
        % taxon_id)

if len(files) == 0:
    print "Didn't find file"
    sys.exit(0)

fh = open(files[0], 'r')
lines = fh.readlines()
fh.close()

accessions = set()
for line in lines:
    if line.strip():
        if line[0] == ">":
            fields = line.split()
            accession = fields[0].split(":")[1]
            accessions.add(accession)

num_covered = 0
for accession in accessions: 
    uniprot_dat_indices = UniProtDatIndex.objects.filter(uniprot_accession = accession) 
    if uniprot_dat_indices:
        uniprot_object = uniprot_dat_indices[0].uniprot
        #if TreeNode.objects.filter(sequence_header__uniprot=uniprot_object, tree__family__active = True).exclude(tree__family__status__exact ='bad'):
        #    num_covered += 1
        if TreeNodeAlignment.objects.filter(sequence_header__uniprot=uniprot_object,tree_node__tree__family__active=True).exclude(tree_node__tree__family__status__exact ='bad'):
            num_covered += 1
        else:
            print "%s is not covered in the database" % accession    
    else:
        print "%s is not even in uniprot_dat_index, so is not covered" % accession
       
taxon = UniProtTaxonomy.objects.get(id = taxon_id)
print "Coverage for %s" % taxon.scientific_name 

print "Source: EBI Reference Proteome (http://www.ebi.ac.uk/reference_proteomes/)"

lineage = taxon.lineage()
print "Taxonomy: %s" % "/".join([str(item) for item in lineage])

print "Number of genes in genome = %d" % len(accessions)

print "Number of genes covered = %d" % num_covered

coverage = float(num_covered) / len(accessions)
print "Coverage = %g" % coverage
all=open ('coverage.csv','a')
all.write ('%s,%s,%s,%s\n' % (taxon.scientific_name,len(accessions),num_covered,coverage))
all.close()
o=open ('%s.coverage' % taxon_id,'w')
o.write ('%s,%s,%s,%s\n' % (taxon.scientific_name,len(accessions),num_covered,coverage))
o.close()
