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
  tree_node = TreeNode.objects.all()[0]
  max_all_left = tree_node._get_max_left_id()
  min_cellular_left = max_all_left
  min_bacterial_left = max_all_left
  min_archaeal_left = max_all_left
  min_eukaryotic_left = max_all_left
  min_noncellular_left = max_all_left
  max_cellular_right = 0
  max_bacterial_right = 0
  max_archaeal_right = 0
  max_eukaryotic_right = 0
  max_noncellular_right = 0
  cellular_left, cellular_right \
      = tree_node._get_left_right_ids('cellular organisms')
  bacterial_left, bacterial_right \
      = tree_node._get_left_right_ids('Bacteria')
  archaeal_left, archaeal_right \
      = tree_node._get_left_right_ids('Archaea')
  eukaryotic_left, eukaryotic_right \
      = tree_node._get_left_right_ids('Eukaryota')
  for header in headers:
    id = header.split()[0]
    id_components = id.split('|')
    if len(id_components) == 3:
      uniprot_accession = id_components[1]
      uniprot = UniProt.objects.get(accession = uniprot_accession)
      left = uniprot.taxon.left_id
      right = uniprot.taxon.right_id
      if left > cellular_left and right < cellular_right:
        # Cellular organism
        if left < min_cellular_left:
          min_cellular_left = left
        if right > max_cellular_right:
          max_cellular_right = right
        # Bacteria?
        if left > bacterial_left and right < bacterial_right:
          if left < min_bacterial_left:
            min_bacterial_left = left
          if right > max_bacterial_right:
            max_bacterial_right = right
        # Archaea?
        elif left > archaeal_left and right < archaeal_right:
          if left < min_archaeal_left:
            min_archaeal_left = left
          if right > max_archaeal_right:
            max_archaeal_right = right
        # Eukaryota?
        elif left > eukaryotic_left and right < eukaryotic_right:
          if left < min_eukaryotic_left:
            min_eukaryotic_left = left
          if right > max_eukaryotic_right:
            max_eukaryotic_right = right
      else:
        # Noncellular organism
        if left < min_noncellular_left:
          min_noncellular_left = left
        if right > max_noncellular_right:
          max_noncellular_right = right

  if min_cellular_left == max_all_left:
    encompassing_taxon = UniProtTaxonomy.objects.filter(
      left_id__lte = min_noncellular_left, 
      right_id__gte = max_noncellular_right).order_by('right_id')[0]
  else:
    # ignore noncellular
    min_kingdom_lefts = [min_bacterial_left, min_archaeal_left,
                        min_eukaryotic_left]
    if len([left for left in min_kingdom_lefts \
            if left <> max_all_left]) > 1:
      # The encompassing taxon would be cellular organisms
      # Return our fake taxon which describes which of the kingdoms
      # is present
      kingdom_strings = []
      if min_bacterial_left <> max_all_left:
        kingdom_strings = kingdom_strings + ['Bacteria']
      if min_archaeal_left <> max_all_left:
        kingdom_strings = kingdom_strings + ['Archaea']
      if min_eukaryotic_left <> max_all_left:
        kingdom_strings = kingdom_strings + ['Eukaryotes']
      scientific_name = ', '.join(kingdom_strings)
      encompassing_taxon \
          = UniProtTaxonomy.objects.get(scientific_name__exact
                                            = scientific_name)
    else:
      encompassing_taxon \
        = UniProtTaxonomy.objects.filter(left_id__lte = min_cellular_left,
            right_id__gte = max_cellular_right).order_by('right_id')[0]

    print encompassing_taxon.__unicode__()
    print ':'.join([t.__unicode__() for t in encompassing_taxon.lineage()])

if __name__ == '__main__':
  main()
