#!/usr/bin/env python

from django.db.models import Q
from pfacts003.phylofacts.models import TreeNode, UniProtTaxonomy, UniProtGO, \
                                        UniProtEC, EC, \
                                        getTreeNodeQuerySetForTreeNodes

from optparse import OptionParser
import sys, re

def GetGenesInTaxonWithOrthologsWithGOAnnotation(uniprot_taxonomy_query,
                                                uniprot_go_query):
  phogs = TreeNode.objects.filter(
            superorthologs__sequence_header__uniprot__go_annotations__in
              = uniprot_go_query)
  canonical_tree_node_of_gene = {}
  phogs_of_gene = {}
  orthologs_of_gene = {}

  for phog in phogs:
    gene_tree_nodes = TreeNode.objects.filter(superorthologous_node = phog,
                  sequence_header__uniprot__taxon__in = uniprot_taxonomy_query)
    if gene_tree_nodes:
      supporting_ortholog_tree_nodes = TreeNode.objects.filter(
                superorthologous_node = phog,
                sequence_header__uniprot__go_annotations__in = uniprot_go_query)
      if supporting_ortholog_tree_nodes:
        for gene_tree_node in gene_tree_nodes:
          gene_id = gene_tree_node.sequence_header.identifier()
          if gene_id not in canonical_tree_node_of_gene:
            canonical_tree_node_of_gene[gene_id] = gene_tree_node
          canonical_tree_node = canonical_tree_node_of_gene[gene_id]
          if canonical_tree_node not in phogs_of_gene:
            phogs_of_gene[canonical_tree_node] = set()
          phogs_of_gene[canonical_tree_node].add(phog)
          if canonical_tree_node not in orthologs_of_gene:
            orthologs_of_gene[canonical_tree_node] = {}
          for supporting_ortholog_tree_node in supporting_ortholog_tree_nodes:
            ortholog_id \
              = supporting_ortholog_tree_node.sequence_header.identifier()
            if ortholog_id not in orthologs_of_gene[canonical_tree_node]:
              orthologs_of_gene[canonical_tree_node][ortholog_id] = set()
            # Store the gene_tree_node, supporting_ortholog_tree_node pair so
            # we can later use them to display the relationship between the
            # gene and its supporting ortholog in a tree, e.g., in Phyloscope
            orthologs_of_gene[canonical_tree_node][ortholog_id].add( 
                      (gene_tree_node, supporting_ortholog_tree_node) )
        
  canonical_tree_nodes \
    = getTreeNodeQuerySetForTreeNodes(canonical_tree_node_of_gene.values())
  return (canonical_tree_nodes, phogs_of_gene, orthologs_of_gene)

def GetGenesInTaxonWithOrthologsWithEC(uniprot_taxonomy_query,
                                                uniprot_ec_query):
  phogs = TreeNode.objects.filter(
            superorthologs__sequence_header__uniprot__ec_annotations__in
              = uniprot_ec_query)
  canonical_tree_node_of_gene = {}
  phogs_of_gene = {}
  orthologs_of_gene = {}

  for phog in phogs:
    gene_tree_nodes = TreeNode.objects.filter(superorthologous_node = phog,
                  sequence_header__uniprot__taxon__in = uniprot_taxonomy_query)
    if gene_tree_nodes:
      supporting_ortholog_tree_nodes = TreeNode.objects.filter(
                superorthologous_node = phog,
                sequence_header__uniprot__ec_annotations__in = uniprot_ec_query)
      if supporting_ortholog_tree_nodes:
        for gene_tree_node in gene_tree_nodes:
          gene_id = gene_tree_node.sequence_header.identifier()
          if gene_id not in canonical_tree_node_of_gene:
            canonical_tree_node_of_gene[gene_id] = gene_tree_node
          canonical_tree_node = canonical_tree_node_of_gene[gene_id]
          if canonical_tree_node not in phogs_of_gene:
            phogs_of_gene[canonical_tree_node] = set()
          phogs_of_gene[canonical_tree_node].add(phog)
          if canonical_tree_node not in orthologs_of_gene:
            orthologs_of_gene[canonical_tree_node] = {}
          for supporting_ortholog_tree_node in supporting_ortholog_tree_nodes:
            ortholog_id \
              = supporting_ortholog_tree_node.sequence_header.identifier()
            if ortholog_id not in orthologs_of_gene[canonical_tree_node]:
              orthologs_of_gene[canonical_tree_node][ortholog_id] = set()
            # Store the gene_tree_node, supporting_ortholog_tree_node pair so
            # we can later use them to display the relationship between the
            # gene and its supporting ortholog in a tree, e.g., in Phyloscope
            orthologs_of_gene[canonical_tree_node][ortholog_id].add( 
                      (gene_tree_node, supporting_ortholog_tree_node) )
        
  canonical_tree_nodes \
    = getTreeNodeQuerySetForTreeNodes(canonical_tree_node_of_gene.values())
  return (canonical_tree_nodes, phogs_of_gene, orthologs_of_gene)

def OutputDBQueryResultsToCSV(canonical_tree_nodes, phogs_of_gene,
                              orthologs_of_gene, csv_file):
  csv_file.write("Gene,Description,PHOGs,Orthologs\n")
  for tree_node in canonical_tree_nodes:
    csv_file.write("%s," % tree_node.sequence_header.identifier())
    csv_file.write('"%s",' % tree_node.sequence_header.description())
    csv_file.write('"%s",' % ','.join([phog.get_accession() for phog in
                                        phogs_of_gene[tree_node]]))
    csv_file.write('"')
    # This is less messy than ','.join()
    first = True
    for ortholog_id in orthologs_of_gene[tree_node].keys():
      if not first:
        csv_file.write(',')
      first = False
      ortholog_tree_node = list(orthologs_of_gene[tree_node][ortholog_id])[0][1]
      csv_file.write("%s (%s)" % (ortholog_id,
        ortholog_tree_node.sequence_header.uniprot.taxon.scientific_name))
    csv_file.write('"')
    csv_file.write('\n')

def main():
  parser = OptionParser(usage='%prog')
  parser.add_option('--taxon', dest='taxon', default='Klebsiella pneumoniae',
                    help="Taxon in which to find genes.")
  parser.add_option('--go_term_type', dest='go_term_type', 
                    default='molecular_function',
                    help="Type of GO term (e.g., molecular_function, " + \
                          "cellular_localization, biological_process).")
  parser.add_option('--go_term', dest='go_term',
    default='low voltage-gated calcium channel activity',
    help="Predicted Gene Ontology (GO) term for which to find genes.")
  parser.add_option('--go_evidence', dest='go_evidence',
    default='EXP',
    help="Evidence code to require for supporting orthologs with the GO term.")
  parser.add_option('--search_type', dest='search_type', default='go',
    help="Type of search: may be 'go' or 'ec', default is 'go'.")
  parser.add_option('--ec', dest='ec', default='2.6.1.82',
  help="Predicted Enzyme Classification (EC) number for which to find genes.")
  parser.add_option('--require_brenda', action='store_true',
    dest='require_brenda', default=True,
    help="Require supporting orthologs with the EC number" + \
          "to be in BRENDA.")
  parser.add_option('--no_require_brenda', action='store_false',
    dest='require_brenda', default=True,
    help="Do not require supporting orthologs with the EC number" + \
          "to be in BRENDA.")
  (options, args) = parser.parse_args()

  tval = options.taxon
  term_type = options.go_term_type
  name = options.go_term
  evcode = options.go_evidence
  search_type = options.search_type.lower()

  # The following code was copied with minor modifications from 
  # pfacts003/dbquery/forms.py

  # Find Taxon
  uptquery = Q(common_name=tval) | \
      Q(scientific_name=tval) | \
      Q(mnemonic=tval)
  try:
      uptquery |= Q(id=int(tval))
  except ValueError:
      pass
  uptquery = UniProtTaxonomy.objects.get(uptquery)
  uptquery = UniProtTaxonomy.objects.filter(
      left_id__gte=uptquery.left_id,
      right_id__lte=uptquery.right_id,
  )
  # End of copied code

  if search_type == 'go':
    # The following code was copied with minor modifications from 
    # pfacts003/dbquery/forms.py
    # Find UniProtGO
    upgquery = UniProtGO.objects.filter(
        go_term__term_type=term_type,
        go_term__name__iexact=name,
    )
    if evcode == 'EXP':
        upgquery = upgquery.filter(go_evidence__id__lte = 6)
    elif evcode != 'ANY':
        upgquery = upgquery.filter(go_evidence__evidence=evcode)
    # End of copied code

    canonical_tree_nodes, phogs_of_gene, orthologs_of_gene \
      = GetGenesInTaxonWithOrthologsWithGOAnnotation(uptquery, upgquery)
  elif search_type == 'ec':
    num_re = re.compile('\d+')
    fields = options.ec.split('.')
    if len(fields) < 4:
      parser.error('Malformed EC number %s' % options.ec)
    try:
      class_number = int(fields[0])
    except ValueError:
      parser.error('Malformed EC number %s' % options.ec)
    ecs = EC.objects.filter(class_number = class_number)
    if fields[1] == '-':
      ecs = ecs.filter(subclass_number__isnull = True)
    else:
      try:
        subclass_number = int(fields[1])
      except ValueError:
        parser.error('Malformed EC number %s' % options.ec)
      ecs = ecs.filter(subclass_number = subclass_number)
      if fields[2] == '-':
        ecs = ecs.filter(subsubclass_number__isnull = True)
      else:
        try:
          subsubclass_number = int(fields[2])
        except ValueError:
          parser.error('Malformed EC number %s' % options.ec)
        ecs = ecs.filter(subsubclass_number = subsubclass_number)
        if fields[3] == '-':
          ecs = ecs.filter(enzyme_number__isnull = True)
        else:
          try:
            enzyme_number = int(fields[3])
            ecs = ecs.filter(enzyme_number = enzyme_number,
                              is_preliminary_f = False)
          except ValueError:
            m = num_re.search(fields[3])
            if m:
              enzyme_number = int(m.group(0))
              ecs = ecs.filter(enzyme_number = enzyme_number,
                                is_preliminary_f = True)
            else:
              parser.error('Malformed EC number %s' % options.ec)
    upequery = UniProtEC.objects.filter(ec__in = ecs)
    if options.require_brenda:
      upequery = upequery.filter(is_in_brenda_f = True)

    canonical_tree_nodes, phogs_of_gene, orthologs_of_gene \
      = GetGenesInTaxonWithOrthologsWithEC(uptquery, upequery)
      
  else:
    parser.error("Unrecognized search type.  Choices are 'go' and 'ec'.")

  OutputDBQueryResultsToCSV(canonical_tree_nodes, phogs_of_gene,
                              orthologs_of_gene, sys.stdout)

if __name__ == '__main__':
  main()
