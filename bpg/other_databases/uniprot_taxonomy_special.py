#!/usr/bin/env python

import sys

import django

from pfacts003.phylofacts import models

def main():
  """Insert PhyloFacts special records into the UniProt taxonomy table"""
  try:
    root_taxon = models.UniProtTaxonomy.objects.get(id__exact = 1)
  except models.UniProtTaxonomy.DoesNotExist:
    # We have not already run on this version of the table
    root_taxon = models.UniProtTaxonomy.objects.create(id=1,
                                                      scientific_name='root')

    # GRJ: What if this changes (i.e., new taxa are found/or there is
    # a re-organization? Couldn't we first check to see if there is a
    # parent, and flag this as an error if a parent already exists?

    # Make most ancestral taxa children of the root,
    # so the NCBI taxonomy will be a tree instead of a forest
    most_ancestral_taxa = models.UniProtTaxonomy.objects.filter(
                      id__in = [131567, # cellular organisms
                                10239,  # Viruses
                                12884,  # Viroids
                                12908,  # unclassified sequences
                                28384,  # other sequences
                              ])
    for taxon in most_ancestral_taxa:
      taxon.parent = root_taxon
      taxon.save()

    # Create taxa (if they don't already exist) representing
    # taxonomic distributions within cellular organisms.
    models.UniProtTaxonomy.objects.get_or_create(id=-3,
                              scientific_name = "Bacteria, Archaea")[0].save()
    models.UniProtTaxonomy.objects.get_or_create(id=-5,
                              scientific_name = "Archaea, Eukaryotes")[0].save()
    models.UniProtTaxonomy.objects.get_or_create(id=-6,
                              scientific_name = "Bacteria, Eukaryotes")[0].save()
    models.UniProtTaxonomy.objects.get_or_create(id=-7,
                              scientific_name = "Bacteria, Archaea, Eukaryotes")[0].save()

    # Give common names to several taxa that have none
    eukaryotes = models.UniProtTaxonomy.objects.get_or_create(
                  scientific_name__exact = "Eukaryota")
    eukaryotes.common_name = "eukaryotes"
    eukaryotes.save()

    slime_molds = models.UniProtTaxonomy.objects.get_or_create(
                  scientific_name__exact = "Mycetozoa")
    slime_molds.common_name = "slime molds"
    slime_molds.save()

    social_amoebae = models.UniProtTaxonomy.objects.get_or_create(
                  scientific_name__exact = "Dictyostelium")
    social_amoebae.common_name = "social amoebae"
    social_amoebae.save()

    euarchontoglires = models.UniProtTaxonomy.objects.get_or_create(
                  scientific_name__exact = "Euarchontoglires")
    euarchontoglires.common_name = \
      "tree shrews, primates, flying lemurs, rodents, and rabbits"
    euarchontoglires.save()
    
if __name__=="__main__":
  main()
