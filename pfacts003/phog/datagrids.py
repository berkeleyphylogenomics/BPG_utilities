"""PHOG DataGrids Module

    The pfacts003.datagrids.datagrids module provides BPG Columns and the
Related Data Grid. This module inherits from RelatedDataGrid to provide the
PHOGDataGrid.
"""

from pfacts003.utils.datagrids.datagrids import *

class OrthologDataGrid(RelatedDataGrid):

  gene = ForeignKeyColumn(field_name="sequence_header", autolink_func = "get_link", sortable=False, 
                label = tooltipped("gene", "Gene ID"))
  uniprot_id = ForeignKeyColumn(field_name="sequence_header.uniprot.uniprot_id",
                      sortable=False, label = tooltipped("uniprot_id", "UniProt ID"))
  species = ForeignKeyColumn(field_name = "sequence_header.uniprot.taxon", autolink_func = "get_link",
              sortable = False, label = tooltipped("species", "Species"))
  description = ForeignKeyColumn(
                  field_name="sequence_header.description",
                  sortable=False, label = tooltipped("description", "Description"))
  in_swissprot = ForeignKeyColumn(field_name="sequence_header.uniprot", autolink_func = "get_swissprot_icon",
                  sortable=False, label = tooltipped("in_swissprot", "Swiss"), expand=False)
  experimental_evidence = ForeignKeyColumn(field_name="sequence_header.uniprot", autolink_func = "get_evidence_icon",
                  sortable=False, label = tooltipped("experimental_evidence", "GO"), expand=False)
  literature = ForeignKeyColumn(field_name="sequence_header.uniprot", autolink_func = "get_literature_icon",
                  sortable=False, label = tooltipped("literature", "Lit."), expand=False)
  family = ForeignKeyColumn(field_name='tree.family', autolink_func = "get_link",
                  sortable=False, 
                  label = tooltipped("family", "PhyloFacts"))
  pfam = PfamColumn(field_name='sequence_header.uniprot', sortable=False, label = tooltipped("pfam", "Pfam<br />Domains"))
  family_type = ForeignKeyColumn(
                  field_name='tree.family.get_alignment_type_str',
                  db_field="family.domain_or_motif_or_full",
                  sortable=False, link = False, label = tooltipped("family_type", "Alignment"))
  phog = PhogColumn(field_name='id', sortable=False, link=False,
                  label = tooltipped("phog", "PhyloFacts<br />Orthology Group"))
# FIXME: RSD 2010/06/15 Commenting out phog_tree and its related
# initializations temporarily, as this requires getting Phyloscope to work.
  phog_tree = PhogTreeColumn(field_name='id',
                  sortable=False, link=True, link_func = link_to_tree,
                  label = tooltipped("phog_tree", "View Tree"))
  ec = ForeignKeyColumn(field_name='sequence_header.uniprot', autolink_func = "get_ec_link",
                  sortable=False, label = tooltipped("ec", "EC"))
# RSD 2010/06/13 Commenting out BioCycColumn for the time being
#  biocyc = BioCycColumn(field_name='sequence_header.uniprot.biocyc_id',
#                  db_field="genbank_uniprot_biocyc.biocyc_id",
#                  sortable=False, expand = False, link = False, label = tooltipped("biocyc", "BioCyc"))
# RSD 2010/06/13 Commenting out this link to the KEGG sequence 
# Will link to the KEGG map instead
#  kegg = KeggColumn(field_name='sequence_header.uniprot.kegg_sequence_id',
#                  db_field="genbank_uniprot_kegg.kegg_id",
#                  sortable=False, expand = False, label = tooltipped("kegg", "KEGG"))
  kegg = ForeignKeyColumn(field_name='sequence_header.uniprot', 
                        autolink_func = "get_kegg_links",
                        sortable = False,
                        label = tooltipped("kegg", "KEGG"))
# RSD 2010/06/12 Commenting out PpiColumn for the time being
#  ppi = PpiColumn(field_name='sequence_header.uniprot', sortable=False, link=False,
#                  label = tooltipped("ppi", "PPI"))

  def __init__(self, request, queryset, extra_context={}, depth=0,
                includePHOGColumns=True, uniprot=None, phog_of_ortholog=None,
                n_sequences_of_family=None,
                ortholog_type = OrthologTypes.SuperOrtholog, threshold=0.0,
                with_family_type=True):
    DataGrid.__init__(self, request, queryset, "")
# RSD 2010/06/13 Commenting out ppi related initializations for the time being
#    self.ppi.uniprot = uniprot
#    self.ppi.ortholog_type = ortholog_type
#    self.ppi.threshold = threshold
    self.phog_tree.uniprot = uniprot
    self.phog.phog_of_ortholog = phog_of_ortholog
    self.phog.ortholog_type = ortholog_type
    self.phog.threshold = threshold
    self.phog_tree.phog_of_ortholog = phog_of_ortholog
    self.phog_tree.n_sequences_of_family = n_sequences_of_family
    self.phog_tree.ortholog_type = ortholog_type
    self.phog_tree.threshold = threshold
    self.listview_template = "phog/ortholog_listview.html"
    self.depth = depth
    self.default_sort = []
    self.extra_context = extra_context
    if includePHOGColumns:
      self.title = \
        mark_safe('Orthologs and In-Species Paralogs from One or More PhyloFacts Families: %s Found' % len(queryset))
      self.default_columns = [
        'gene',
        'species',
        'description',
        'in_swissprot', 
        'experimental_evidence',
        'ec',
        'kegg',
        'biocyc',
        'ppi',
        'pfam',
        'family', 
        ]+ (with_family_type and ['family_type'] or []) +[
        'phog',
        'phog_tree',
        'literature',
        ]
    else:
      self.title = 'PhyloFacts Orthology Group Members'
      self.default_columns = [
        'gene',
        'species',
        'description',
        'in_swissprot', 
        'experimental_evidence',
        'ec',
        'kegg',
        'biocyc',
        'ppi',
        'literature',
        ]

class ApproximateMatchesDataGrid(RelatedDataGrid):
  description = ForeignKeyColumn(
                  field_name="sequence_header.description",
                  sortable=False, label = tooltipped("description", "Description"))
  species = ForeignKeyColumn(field_name = "sequence_header.uniprot.taxon", autolink_func = "get_link",
              sortable = False, label = tooltipped("species", "Species"))
  gene = ForeignKeyColumn(field_name="sequence_header", autolink_func = "get_link", sortable=False, 
                label = tooltipped("gene", "Gene ID"))
  percent_id = PercentIdColumn(field_name="sequence_header.id", sortable=False, 
                      css_class = 'percent_id',
                      label = tooltipped("percent_id", "%ID"))
  percent_cover = PercentCoverColumn(field_name="sequence_header.id", sortable=False, 
                      css_class = 'percent_cover',
                      label = tooltipped("percent_cover", "%Cover"))
  orthologs = OrthologsColumn(field_name='sequence_header',
                  sortable=False, label = tooltipped("orthologs", "View Orthologs"))
  alignment_length = AlignmentLengthColumn(field_name="sequence_header.id", sortable=False, 
                      css_class = 'alignment_length',
                      label = tooltipped("aligned_length", "Align.<br />Len."))
  e_value = EValueColumn(field_name="sequence_header.id", sortable=False, 
                      css_class = 'e_value',
                      label = tooltipped("e_value", "BLAST<br />E-Value"))
  blast_score = BlastScoreColumn(field_name="sequence_header.id", sortable=False, 
                      css_class = 'blast_score',
                      label = tooltipped("blast_score", "BLAST<br />Score"))
  blast_alignment = BlastAlignmentColumn(field_name="sequence_header.id", sortable=False, 
                      label = tooltipped("blast_alignment", "View<br />Alignment"))
  def precompute_objects(self):
    RelatedDataGrid.precompute_objects(self)

    # Sort QuerySet objects by BLAST score in dictionary
    # This scheme will only work if all objects are on a single page
    objects = []    
    for obj in self.page.object_list:
      key = str(obj.sequence_header.id)
      objects.append((self.protseq_alignment_dict[key]['score'], obj)) 
    objects.sort()
    objects.reverse()
    self.page.object_list = [item[1] for item in objects] 

    # Now reconstruct rows in the table from the sorted objects
    self.rows = []
    for obj in self.page.object_list:
      self.rows.append({
        'object': obj,
        'cells': [column.render_cell(obj) for column in self.columns]
      })
      
  def __init__(self, request, queryset, protseq_alignment_dict, extra_context={}, depth=0):
    DataGrid.__init__(self, request, queryset, "BLAST Matches")
    self.listview_template = "phog/ortholog_listview.html"
    self.paginate_by = 100
    self.depth = depth
    self.default_sort = []
    self.extra_context = extra_context
    self.protseq_alignment_dict = protseq_alignment_dict
    self.default_columns = [
      'gene',
      'orthologs',
      'species',
      'description',
      'percent_id',
      'e_value',
      'blast_score',
      ]
    try:
      self.blast_results_file_name = self.extra_context['blast_results_file_name']
      self.default_columns.append('blast_alignment')
    except KeyError, e:
      pass
#class PHOG2PHOGDataGrid(RelatedDataGrid):
#  species = ForeignKeyColumn(field_name = "get_taxon", autolink_func = "get_link",
#              sortable = False, label = tooltipped("species", "Species"))
#  gene1 = InteractingGeneColumn(field_name = "id", sortable = False,
#                  label = tooltipped("gene 1", "Gene 1"))
#  gene2 = InteractingGeneColumn(field_name = "id", sortable = False,
#                  label = tooltipped("gene 2", "Gene 2"))
#  source = Column(field_name = "id", sortable = False, link = True,
#                  link_func = DataGrid.link_to_object,
#                  label = tooltipped("source", "Source"))
#  def __init__(self, request, queryset, uniprot1_of_edge, uniprot2_of_edge,
#              extra_context= {}, depth=0):
#    DataGrid.__init__(self, request, queryset, "Observed Interactions")
#    self.listview_template = "servers/phog/ortholog_listview.html"
#    self.depth = depth
#    self.default_sort = []
#    self.extra_context = extra_context
#    self.gene1.uniprot_of_edge = uniprot1_of_edge
#    self.gene2.uniprot_of_edge = uniprot2_of_edge
#    self.default_columns = [
#      'species',
#      'gene1',
#      'gene2',
#      'source',
#    ] 
