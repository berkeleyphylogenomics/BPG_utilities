#!/usr/bin/python

from Bio.SeqUtils import CheckSum
import urllib2
import glob, sys, os.path
from pfacts003.phylofacts.models import *
import logging

'''
This script accepts a taxon id, finds the file with the genome for that
taxon in the /clusterfs/ohana/external/genomes/QuestForOrthologs/Release5/
directory, figures out if each accession is found in any family, and
prints out  a report on the coverage in PhyloFacts for the genome.
'''
def setup_logger(taxon_id):
  # Set logging so that replacing print statements with self.log
  # writes logs to a file
  logger = logging.getLogger()
  hdlr = logging.FileHandler('%s.log' % taxon_id)
  formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
  hdlr.setFormatter(formatter)
  logger.addHandler(hdlr)
  logger.setLevel(logging.INFO)
  return  logger

def log(logger,string):
  '''Actual print statement replacement.'''
  logger.info(string)

def get_family_size(TNA):
  '''get the size of family, and count how many of them has over 3 sequences'''
  max_size=0
  families=set([tree_node_alignment.tree_node.tree.family \
            for tree_node_alignment in TNA])
  for family in families: 
    size =len (family.canonical_root_node().get_included_leaves())
    if max_size <size:
      max_size=size

  if max_size >=3:
    return 1
  else:
    return 0

  
def main():
  if  len(sys.argv) < 2:
      print "usage: specify taxonomy id"
      sys.exit(1)

  taxon_id = sys.argv[1]
  input_file = sys.argv[2]
  
  logging.basicConfig()
  logger=setup_logger(taxon_id)
  files = glob.glob("/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5/%s_*.fasta" 
          % taxon_id)

  if len(files) == 0:
      files=glob.glob("/clusterfs/ohana/bpg/coverage/redundant/pfam/after_17GHG/ID/QFO/%s_*.fasta" % taxon_id)
      if len(files) == 0:
          #print "Didn't find file"
          fh = open (input_file, 'r')
          #sys.exit(0)
  
  #fh = open(files[0], 'r')
  lines = fh.readlines()
  fh.close()

  accessions = set()
  for line in lines:
      if line.strip():
          if line[0] == ">":
              fields = line.split()
              try:
                accession = fields[0].split(":")[1]
              except:
                accession = fields[0].split("|")[1]
              accessions.add(accession)

  o=open ('%s.coverage' % taxon_id,'w')
  num_covered = GHG_num_covered = GHG_over2_covered = Pfam_num_covered=0
  for accession in accessions: 
      # First search the uniprot accession in the UniProt_Dat_Index table
      uniprot_dat_indices = UniProtDatIndex.objects.filter(uniprot_accession = accession) 
      if uniprot_dat_indices:
          uniprot_object = uniprot_dat_indices[0].uniprot
          if TreeNodeAlignment.objects.filter(sequence_header__uniprot=uniprot_object,tree_node__tree__family__active=True).exclude(tree_node__tree__family__status__exact='bad'):
            num_covered +=1
            TNA=TreeNodeAlignment.objects.filter(sequence_header__uniprot=uniprot_object,tree_node__tree__family__active=True,tree_node__tree__family__family_type='G').exclude(tree_node__tree__family__status__exact='bad')
            if len(TNA)>0:
                GHG_num_covered += 1
                GHG_over2_covered += get_family_size(TNA)
            if TreeNodeAlignment.objects.filter(sequence_header__uniprot=uniprot_object,tree_node__tree__family__active=True,tree_node__tree__family__family_type='C').exclude(tree_node__tree__family__status__exact='bad'):
                Pfam_num_covered += 1
          else:
               log(logger,"%s is not covered in the database\n" % accession)    
      # if the accession is not in the uniprot_dat_index table, use seguid to
      # find identical sequences
      else:
          log(logger,"%s is not in uniprot_dat_index, try the seguid\n" % accession)
          uniprot_accession = accession
          if uniprot_accession_re1.match(uniprot_accession) or \
            uniprot_accession_re2.match(uniprot_accession):
              fasta_file='%s.fasta' % uniprot_accession
              cmd = 'wget http://www.uniprot.org/uniprot/%s' % fasta_file 
              try:
                  os.system(cmd)
              except:
                  log(logger,"Unable to download sequence from UniProt\n")
              response=open(fasta_file,'r')
              record = SeqIO.parse(response, 'fasta').next()
              response.close()
              os.remove(fasta_file)
              seguid = CheckSum.seguid(record.seq) 
              sequence_objects = Sequence.objects.filter(seguid__exact = seguid)
              if sequence_objects:
                  if TreeNodeAlignment.objects.filter(sequence_header__sequence__in = sequence_objects,tree_node__tree__family__active=True).exclude(tree_node__tree__family__status__exact='bad'):
                      num_covered += 1
                      TNA=TreeNodeAlignment.objects.filter(sequence_header__sequence__in=sequence_objects,tree_node__tree__family__active=True,tree_node__tree__family__family_type='G').exclude(tree_node__tree__family__status__exact='bad')
                      if len(TNA)>0:
                          GHG_num_covered += 1
                          GHG_over2_covered += get_family_size(TNA)
                      if TreeNodeAlignment.objects.filter(sequence_header__sequence__in=sequence_objects,tree_node__tree__family__active=True,tree_node__tree__family__family_type='C').exclude(tree_node__tree__family__status__exact='bad'):
                          Pfam_num_covered += 1
                  else:
                      o.write("There are no families containing %s.\n" % accession)
              else:
                   log(logger,"%s is not in the PhyloFacts 3 database.\n" % accession)
          else:
              print "The argument must be a valid UniProt accession\n" 
                  
  taxon = UniProtTaxonomy.objects.get(id = taxon_id)
  print "Coverage for %s" % taxon.scientific_name 

  print "Source: EBI Reference Proteome (http://www.ebi.ac.uk/reference_proteomes/)"


  lineage = taxon.lineage()
  print "Taxonomy: %s" % "/".join([str(item) for item in lineage])

  print "Number of genes in genome = %d" % len(accessions)

  print "Number of genes covered = %d" % num_covered

  coverage_GHG = float(GHG_num_covered) / len(accessions)
  coverage_GHG_over2 = float(GHG_over2_covered) / len(accessions)
  coverage_Pfam = float(Pfam_num_covered) / len(accessions)
  coverage_any = float(num_covered) / len(accessions)

  log(logger,"Coverage = %g" % coverage_any)
  o.write ('taxon ID,scientific name,# of proteins,# covered,% covered,# covered by \
  GHG,% covered by GHG,# covered by GHG of size>=3,% covered by GHG of size >=3,# covered by Pfam, %covered by Pfam\n')
  o.write ('%s,%s,%d,%d,%3.1f,%d,%3.1f,%d,%3.1f,%d,%3.1f\n'
%(taxon_id,taxon.scientific_name,len(accessions),num_covered,coverage_any*100,GHG_num_covered,coverage_GHG*100,GHG_over2_covered,coverage_GHG_over2*100,Pfam_num_covered,coverage_Pfam*100))
  o.close()

if __name__=='__main__':
  main()
