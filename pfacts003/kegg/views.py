from django.shortcuts import render_to_response
from pfacts003.kegg.forms import *
from pfacts003.phylofacts.models import OrthologTypes, KEGG_Map, KEGG_Map_EC, \
EC, UniProtEC, UniProtTaxonomy, TreeNode, getTreeNodeQuerySetForTreeNodes
import string
from django.db.models import Q

ec_translation = string.maketrans('.-','__')

def kegg_map(request, kegg_map_id, taxon_id=None):
  context = {}
  kegg_map = KEGG_Map.objects.get(id = int(kegg_map_id[3:]))
  context['kegg_map_id'] = kegg_map_id
  context['kegg_map_number'] = kegg_map_id[3:]
  context['kegg_map_title'] = kegg_map.title
  if taxon_id != None:
    newGET = request.GET.copy()
    newGET.update({'taxon_id': taxon_id})
    request.GET = newGET
    context['taxon_id'] = taxon_id
    form = KeggForm(request.GET)
    form.taxon_id = taxon_id
  elif 'taxon_id' in request.GET:
    taxon_id = request.GET['taxon_id'].rstrip()
    context['taxon_id'] = request.GET['taxon_id']
    form = KeggMapForm(request.GET)
    form.taxon_id = request.GET['taxon_id']
  else:
    form = KeggMapForm()
  context['form'] = form
  if form.is_bound and form.is_valid() and taxon_id != None and \
      taxon_id != u'':
    if form.non_field_errors():
      context['form_error'] = form.non_field_errors().as_text()
    else:
      # For some reason, ortholog_type comes in as Unicode
      try:
        ortholog_type = int(form.cleaned_data['ortholog_type']) 
      except ValueError:
        ortholog_type = OrthologTypes.SuperOrtholog
      threshold = form.cleaned_data['threshold']
      kegg_map_ecs = KEGG_Map_EC.objects.filter(kegg_map = kegg_map)
      ecs = [kegg_map_ec.ec for kegg_map_ec in kegg_map_ecs]
      require_brenda = bool(form.cleaned_data['require_brenda'])
      taxon = UniProtTaxonomy.objects.get(id = taxon_id)
      nodes_to_color = {}
      phogs_of_uniprot_of_ec = {}
      for ec in ecs:
        if require_brenda:
          uniprot_ecs = UniProtEC.objects.filter(ec = ec,
                                                is_in_brenda_f = True)
        else:
          uniprot_ecs = UniProtEC.objects.filter(ec = ec)
        uniprots = [uniprot_ec.uniprot for uniprot_ec in uniprot_ecs]
        query_nodes = TreeNode.objects.filter(
                sequence_header__uniprot__in=uniprots)
        phogs = set()
        for query_node in query_nodes:
            containing_phog \
              = query_node.get_containing_phog(threshold=threshold)
            if containing_phog:
              phogs.add(containing_phog)
        phogs_of_uniprot_of_ec[ec] = {}
        for phog in phogs:
          leaves = phog.get_included_leaves(ortholog_type = ortholog_type,
                                              threshold = threshold)
          orthologs = leaves.exclude(
                sequence_header__uniprot__taxon__left_id__lt = taxon.left_id
                ).exclude(
              sequence_header__uniprot__taxon__right_id__gt = taxon.right_id)
          known_leaves = leaves.filter(sequence_header__uniprot__in = uniprots)
          for ortholog in orthologs:
            if ortholog.sequence_header.uniprot \
                not in phogs_of_uniprot_of_ec[ec]:
              phogs_of_uniprot_of_ec[ec][ortholog.sequence_header.uniprot] = {}
              phogs_of_uniprot_of_ec[ec][ortholog.sequence_header.uniprot][
                'phog_accessions'] = set()
              phogs_of_uniprot_of_ec[ec][ortholog.sequence_header.uniprot][
                                    'known_uniprots'] = set()
            phogs_of_uniprot_of_ec[ec][ortholog.sequence_header.uniprot][
                                    'phog_accessions'].add(phog.get_accession(
                                    ortholog_type = ortholog_type, 
                                    threshold = threshold))
            phogs_of_uniprot_of_ec[ec][ortholog.sequence_header.uniprot][
                                    'known_uniprots'] \
            = phogs_of_uniprot_of_ec[ec][ortholog.sequence_header.uniprot][
                                    'known_uniprots'] \
              | set([leaf.sequence_header.uniprot for leaf in known_leaves])
        if len(phogs_of_uniprot_of_ec[ec]) > 0:
          nodes_to_color[ec.__str__().replace('.','_').replace('-','_')] = True
      context['phogs_of_uniprot_of_ec'] = phogs_of_uniprot_of_ec
      context['nodes_to_color'] = nodes_to_color
  return render_to_response('kegg/%s.html' % kegg_map_id, context)
