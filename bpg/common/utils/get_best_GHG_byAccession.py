#!/usr/bin/python
from Bio.SeqUtils import CheckSum
import urllib2
import glob, sys, os.path
from pfacts003.phylofacts.models import *

def get_all_taxon_ids(uniq_families): # get the taxonomic ID from each sequence within a family
    max_taxon_num=0
    max_family=''
    for family in uniq_families:
        family_score = family.get_score()
        if (family_score > max_taxon_num):
            max_family = family
            max_taxon_num = family_score
    return max_family,max_taxon_num

def main():

# This script takes UniProt Accession as input and finds the best GHG books
# that contains the Accession. Best book means having the most diverse taxa.
# Note: In PhyloFacts3 there are families without trees (normally because all
# the sequences in the family are identical, which makes the tree has no branch
# length), this script cannot find these families.
    if  len(sys.argv) < 2:
        print "usage: %prog UniProt_accession outputfile "
        sys.exit(1)

    accession = sys.argv[1]
# find the UniProt index for the accession
    uniprot_dat_indices = UniProtDatIndex.objects.filter(uniprot_accession = accession)
    if uniprot_dat_indices:
        uniprot_object = uniprot_dat_indices[0].uniprot
        TreeNodeAlignments=TreeNodeAlignment.objects.filter(sequence_header__uniprot=uniprot_object,tree_node__tree__family__active=True,tree_node__tree__family__family_type='G').exclude(tree_node__tree__family__status__exact ='bad')
        if not len(TreeNodeAlignments)==0:
            families=set([tree_node_alignment.tree_node.tree.family \
              for tree_node_alignment in TreeNodeAlignments])
            (max_family_id, max_taxa_num)=get_all_taxon_ids(families)
            print ("%s,%s,%s" %(accession,max_family_id,max_taxa_num))
        else:
            print ("%s is not covered in the database\n" % accession)
    else:
# find the seguid from UniProt
        print ("%s is not in uniprot_dat_index, try the seguid\n" % accession)
        uniprot_accession = accession
        if uniprot_accession_re1.match(uniprot_accession) or \
          uniprot_accession_re2.match(uniprot_accession):
            try:
                response = urllib2.urlopen('http://www.uniprot.org/uniprot/%s.fasta'
                                % uniprot_accession)
            except urllib2.HTTPError:
                print ("Unable to download sequence from UniProt\n")
            record = SeqIO.parse(response, 'fasta').next()
            seguid = CheckSum.seguid(record.seq)
            sequence_objects = Sequence.objects.filter(seguid__exact = seguid)
            if sequence_objects:
                TreeNodeAlignments=TreeNodeAlignment.objects.filter \
                 (sequence_header__sequence__in = sequence_objects, \
                 tree_node__tree__family__active=True, \
                 tree_node__tree__family__family_type='G') \
                 .exclude(tree_node__tree__family__status__exact='bad')
                if not len(TreeNodeAlignments)==0:
                    families=set([tree_node_alignment.tree_node.tree.family \
  for tree_node_alignment in TreeNodeAlignments])
                    (max_family_id,max_taxa_num)=get_all_taxon_ids(families)
                    print ("%s,%s,%s" % (accession,max_family_id,max_taxa_num))
                else:
                    print ("There are no families containing %s.\n" % accession)
            else:
                print ("%s is not in the PhyloFacts 3 database.\n" % accession)
        else:
            print ("The argument must be a valid UniProt accession\n")

if __name__ == '__main__':
    main()
