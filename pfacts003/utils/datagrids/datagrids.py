"""BPG Columns and Datagrids

    This module provides BPG specific columns and the RelatedDataGrid. You
should not add your own DataGrids to this file. Instead, you create a separate
file, import the namespaces from this file, and subclass RelatedDataGrid. 

from pfacts003.utils.datagrids.datagrids import *

class MyNewDataGrid(RelatedDataGrid):
    id = PfamColumn("id", sortable=False)
    ...
"""

import string

from django.core.paginator import InvalidPage, QuerySetPaginator
from django.conf import settings
from django.http import Http404
from django.template.context import RequestContext
from django.template.loader import render_to_string
from django.utils.safestring import mark_safe
from djblets.datagrid.grids import Column, DataGrid

from pfacts003.utils.links import make_link, make_tooltipped_text, \
                  make_pfam_links, make_uniprot_pfam_links, \
                  alignment_blast_pairwise_link, biocyc_icon, kegg_icon, \
                  dip_link, uniprot_lit_icon, uniprot_orthologs_link, \
                  uniprot_link, family_orthoscope_link_in_phog_table, \
                  phog_tree_link
from pfacts003.phylofacts.models import Family, TreeNode, OrthologTypes

def link_to_tree(obj, value):
  return value.get_tree_url()


##############
# Tooltips
##############

# Temporarily removed because it has, of yet, not been needed.
# If any of the following Columns reference tooltips, this should
# be uncommented.
## Tooltips are selected by a tag 'tip tag' defined in
## /media/js/phog_tooltip.js.
# Tip tags that begin with an underscore are used for the column headings.

def tooltipped(column, heading):
  return make_tooltipped_text(heading, '_'+column)


#########
# Columns
#########

# Defined first for dependencies
class ForeignKeyColumn(Column):

  def __init__(self, autolink_func = '', *args, **kwargs):
    Column.__init__(self, *args, **kwargs)
    self.autolink_func = autolink_func

  def render_data(self, obj):
    """
    Renders the column data to a string.
    Expands out multiple dots in the field name as foreign keys.
    """
    value = obj
    for field_component in self.field_name.split('.'):
      value = getattr(value, field_component)
      if value is None:
        return 'N/A'

    if callable(value):
      return value()
    else:
      if self.autolink_func:
        try:
          func = getattr(value, self.autolink_func)
        except AttributeError:
          return 'Error &ndash; Object &ldquo;%s&rdquo; lacks rendering function &ldquo;%s&rdquo;' \
              % (value.__class__.__name__, self.autolink_func)
        value = func()
      return value

class AlignmentLengthColumn(ForeignKeyColumn):
  def render_data(self, obj):
    protseq_id = ForeignKeyColumn.render_data(self, obj)
    return self.datagrid.protseq_alignment_dict[str(protseq_id)]['alignment_length']

class BioCycColumn(ForeignKeyColumn):
  def render_data(self, obj):
    return biocyc_icon(ForeignKeyColumn.render_data(self, obj))

class BlastAlignmentColumn(ForeignKeyColumn):
  def render_data(self, obj):
    """
    Renders the column data as "Alignment" with link to BLAST pairwise
    alignment page
    """
    value = "bpgprotseq%09d" % ForeignKeyColumn.render_data(self, obj)
    return alignment_blast_pairwise_link('Alignment', value, self.datagrid.blast_results_file_name)


class BlastScoreColumn(ForeignKeyColumn):
  def render_data(self, obj):
    protseq_id = ForeignKeyColumn.render_data(self, obj)
    return '%d' % self.datagrid.protseq_alignment_dict[str(protseq_id)]['score']


class EValueColumn(ForeignKeyColumn):
  def render_data(self, obj):
    protseq_id = ForeignKeyColumn.render_data(self, obj)
    e_value = self.datagrid.protseq_alignment_dict[str(protseq_id)]['e_value']
    if e_value == 0:
      return e_value
    else:
      return '%.1e' % e_value

class InteractingGeneColumn(ForeignKeyColumn):
  def __init__(self, *args, **kwargs):
    ForeignKeyColumn.__init__(self, *args, **kwargs)
    self.uniprot_of_edge = None
  def render_data(self, obj):
    edge = ForeignKeyColumn.render_data(self, obj)
    if self.uniprot_of_edge and edge in self.uniprot_of_edge:
      uniprot = self.uniprot_of_edge[edge]
      return uniprot_orthologs_link(uniprot.uniprot_id, text = uniprot.uniprot_id)

class InteractionSourceColumn(ForeignKeyColumn):
  def __init__(self, *args, **kwargs):
    ForeignKeyColumn.__init__(self, *args, **kwargs)
    self.interaction_queryset = None
  def render_data(self, obj):
    own_uniprot_id = ForeignKeyColumn.render_data(self, obj)
    if own_uniprot_id == '' or not self.interaction_queryset:
      return ''
    dip_interactions = self.interaction_queryset.filter(
      dip_nodes__dip_node__uniprots__uniprot__uniprot_id__exact = own_uniprot_id)
    return dip_link(dip_interactions[0].id)


class KeggColumn(ForeignKeyColumn):
  def render_data(self, obj):
    return kegg_icon(ForeignKeyColumn.render_data(self, obj))


class OrthologsColumn(Column):
  def render_data(self, obj):
    return Column.render_data(self, obj).get_orthologs_link()


class PercentCoverColumn(ForeignKeyColumn):
  def render_data(self, obj):
    protseq_id = ForeignKeyColumn.render_data(self, obj)
    return '%.1f' % self.datagrid.protseq_alignment_dict[str(protseq_id)]['percent_coverage']

class PercentIdColumn(ForeignKeyColumn):
  def render_data(self, obj):
    protseq_id = ForeignKeyColumn.render_data(self, obj)
    return '%.1f' % self.datagrid.protseq_alignment_dict[str(protseq_id)]['percent_identity']

class PfamColumn(ForeignKeyColumn):
  def render_data(self, obj):
    uniprot = ForeignKeyColumn.render_data(self, obj)
    if uniprot == 'N/A':
        return 'N/A'
    return make_uniprot_pfam_links(uniprot)

class PhogColumn(ForeignKeyColumn):
  def __init__(self, *args, **kwargs):
    ForeignKeyColumn.__init__(self, *args, **kwargs)
    self.phog_of_ortholog = None
  def render_data(self, obj):
    """
    Renders the column data as a link to the PHOG page
    """
    own_node = ForeignKeyColumn.render_data(self, obj)
    if not own_node:
      return ''
    own_node = TreeNode.objects.get(id__exact = own_node)
    if self.phog_of_ortholog:
      if own_node in self.phog_of_ortholog:
        phog_node = self.phog_of_ortholog[own_node]
      else:
        phog_node = own_node.superorthologous_node
    else:
      phog_node = own_node.superorthologous_node
    if phog_node:
      url = phog_node.get_absolute_url(self.ortholog_type, self.threshold)
      link = make_link(phog_node.get_accession(self.ortholog_type, self.threshold), 
                       url, 'phog')
      return link
    else:
      return ''

class PhogTreeColumn(ForeignKeyColumn):
  def __init__(self, *args, **kwargs):
    ForeignKeyColumn.__init__(self, *args, **kwargs)
    self.uniprot = None
    self.phog_of_ortholog = None
    self.ortholog_type = OrthologTypes.SuperOrtholog
    self.threshold = 0.0
    self.n_sequences_of_family = None
  def render_data(self, obj):
    """
    Renders the column data as "Full PHOG" with links to the Full and PHOG trees
    """
    own_node = ForeignKeyColumn.render_data(self, obj)
    if not own_node:
      return ''
    own_node = TreeNode.objects.get(id__exact = own_node)
    if self.phog_of_ortholog:
      if own_node in self.phog_of_ortholog:
        phog_node = self.phog_of_ortholog[own_node]
      else:
        phog_node = own_node.superorthologous_node
    else:
      phog_node = own_node.superorthologous_node
    if phog_node:
      family = phog_node.tree.family.__unicode__();
      n_sequences = self.n_sequences_of_family[phog_node.tree.family]
      ortholog_leftid_str = str(obj.left_id)
      query_leftid_str = ""
      if self.uniprot:
        query_node_set = TreeNode.objects.filter(tree__exact = phog_node.tree,
                                  sequence_header__uniprot__exact = self.uniprot)
        if query_node_set:
          query_node = query_node_set[0]
          query_leftid_str = str(query_node.left_id)
      link = "%s&nbsp;%s" % (
          family_orthoscope_link_in_phog_table(family, n_sequences),
          phog_tree_link(phog_node.get_accession(self.ortholog_type, self.threshold), 
            ortholog_leftid_str, query_leftid_str))
      return link
    else:
      return ''

class UniProtColumn(Column):
  def render_data(self, obj):
    value = Column.render_data(self, obj)
    return uniprot_link(value)

class UniProtOrthologsColumn(ForeignKeyColumn):
  def render_data(self, obj):
    value = ForeignKeyColumn.render_data(self, obj)
    return uniprot_orthologs_link(value)

  def render_listview(self):
      """
      Renders the standard list view of the grid.

      This can be called from templates.

      Just like the method in the superclass, except it updates the 
      context with self.extra_context first.
      """
      self.load_state()

      context = RequestContext(self.request, {
                  'datagrid': self,
                  'is_paginated': self.page.has_other_pages(),
                  'results_per_page': self.paginate_by,
                  'has_next': self.page.has_next(),
                  'has_previous': self.page.has_previous(),
                  'page': self.page.number,
                  'next': self.page.next_page_number(),
                  'previous': self.page.previous_page_number(),
                  'last_on_page': self.page.end_index(),
                  'first_on_page': self.page.start_index(),
                  'pages': self.paginator.num_pages,
                  'extraquery': 'search_key=blah',
                  'hits': self.paginator.count,
                  'page_range': self.paginator.page_range,
                })
      context.update(self.extra_context)
      return mark_safe(render_to_string(self.listview_template, context))


#########################################
# Data Grids
#########################################

class RelatedDataGrid(DataGrid):
  """Version of djblets DataGrid modified for BPG needs

  This class is almost identical to it's superclass, DataGrid. The only
differences are that the precomputer_objects and render_listview methods have
been overridden. However, the code in both methods is very similar to the
DataGrid superclass code for these methods. Small, customized 'tweaks' have
been made to the original and are documented per method.
  """

  def precompute_objects(self):
      """
      Builds the queryset and stores the list of objects for use in
      rendering the datagrid.

      Just like the method in the superclass, 
      except we use self.depth as the depth for select_related
      """
      query = self.queryset
      use_select_related = False

      # Generate the actual list of fields we'll be sorting by
      sort_list = []
      for sort_item in self.sort_list:
          if sort_item[0] == "-":
              base_sort_item = sort_item[1:]
              prefix = "-"
          else:
              base_sort_item = sort_item
              prefix = ""

          if sort_item and base_sort_item in self.db_field_map:
              db_field = self.db_field_map[base_sort_item]
              sort_list.append(prefix + db_field)

              # Lookups spanning tables require that we query from those
              # tables. In order to keep things simple, we'll just use
              # select_related so that we don't have to figure out the
              # table relationships. We only do this if we have a lookup
              # spanning tables.
              if '.' in db_field:
                  use_select_related = True

      if sort_list:
          query = query.order_by(*sort_list)

      if use_select_related:
          query = query.select_related(depth=self.depth)

      self.paginator = QuerySetPaginator(query, self.paginate_by,
                                       self.paginate_orphans)

      page_num = self.request.GET.get('page', 1)

      # Accept either "last" or a valid page number.
      if page_num == "last":
          page_num = self.paginator.num_pages

      try:
          self.page = self.paginator.page(page_num)
      except InvalidPage:
          raise Http404

      self.rows = []

      for obj in self.page.object_list:
          self.rows.append({
              'object': obj,
              'cells': [column.render_cell(obj) for column in self.columns]
          })

  def render_listview(self):
      """
      Renders the standard list view of the grid.

      This can be called from templates.

      Just like the method in the superclass, except it updates the 
      context with self.extra_context first.
      """
      self.load_state()

      context = RequestContext(self.request, {
                  'datagrid': self,
                  'is_paginated': self.page.has_other_pages(),
                  'results_per_page': self.paginate_by,
                  'has_next': self.page.has_next(),
                  'has_previous': self.page.has_previous(),
                  'page': self.page.number,
                  'next': self.page.next_page_number(),
                  'previous': self.page.previous_page_number(),
                  'last_on_page': self.page.end_index(),
                  'first_on_page': self.page.start_index(),
                  'pages': self.paginator.num_pages,
                  'hits': self.paginator.count,
                  'page_range': self.paginator.page_range,
                })
      context.update(self.extra_context)
      return mark_safe(render_to_string(self.listview_template, context))
