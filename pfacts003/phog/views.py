import csv, string, glob, os
import xml.sax.saxutils as saxutils

from django.template import Context, RequestContext, loader
from django.shortcuts import render_to_response
from django.http import QueryDict, HttpResponse, HttpResponseRedirect
from django.forms.util import ErrorList
from django.utils.safestring import mark_safe

from pfacts003.phog.forms import *
from pfacts003.phog.datagrids import OrthologDataGrid, \
     ApproximateMatchesDataGrid 
from pfacts003.phog.orthologs import OrthologTypes, getOrthologQuerySet, \
      getPHOGQuerySet
from pfacts003.phylofacts.models import TreeNode, OrthologTypes, \
    get_ortholog_type_threshold_from_phog_accession, UniProtTaxonomy, Family
from pfacts003.utils.get_identifier_from_fasta import \
                 get_uniprot_id_from_fasta, \
                 get_approximately_matching_sequence_query_set, \
                 get_alignment_blast_pairwise
from pfacts003.utils.json_data_objects import annotatedTree, \
          superorthologousNodes, \
          hyperNeighborhood, simpleHyperNeighborhood, getLineage, getMRCA
from pfacts003.utils.messages import not_in_phog
from pfacts003.utils.links import make_pfam_links, make_link
from pfacts003.phylofacts.models import UniProt, TreeNodeAlignment
from pfacts003.phylofacts.models import PHOGRow

trivial_translation = string.maketrans('', '')
dotlowercase = '.' + string.lowercase

def phog_url(phog, ortholog_type = OrthologTypes.SuperOrtholog,
                  threshold = 0.0):
  return phog.get_absolute_url(ortholog_type, threshold)


def phog_link(phog, ortholog_type = OrthologTypes.SuperOrtholog,
                    threshold = 0.0):
  return make_link(phog.get_accession(ortholog_type, threshold),
                  phog_url(phog, ortholog_type, threshold), 'phog')


def get_parameter_from_GET(parameter, request):
  if parameter in request.GET:
    return request.GET[parameter]
  else: 
    return ""

def index(request):
    t = loader.get_template('base_templates/index.html')
    c = Context({})
    return HttpResponse(t.render(c))

def ajax_get_lineage(request):
    return HttpResponse(getLineage(request.GET.get('taxid')),
        mimetype = 'application/javascript') # Avoid security problems

def ajax_get_mrca(request):
    return HttpResponse(getMRCA(request.GET.getlist('taxid')), # Note: getlist()
        mimetype = 'application/javascript') # Avoid security problems

def alignment_blast_pairwise(request):
    t = loader.get_template('phog/alignment_blast_pairwise.html')
    c = RequestContext(request, {
        'identifier': get_parameter_from_GET('identifier', request),
        'location': get_parameter_from_GET('location', request),
        'file': get_parameter_from_GET('file', request)})
    get_alignment_blast_pairwise(c)
    return HttpResponse(t.render(c))
  
def front(request):
    t = loader.get_template('main/front.html')
    c = RequestContext(request, {})
    return HttpResponse(t.render(c))

def file_listing(request):
    family_name = get_parameter_from_GET('family', request)
    if False and not 'file' in request.GET:
      dir = family_dir(family_name, webdir = True)
    else:
      file_name = family_nj_file_name(family_name, webdir = True)
      content = ''
      try:
        file = open(file_name, 'r')
        line_num = 0
        content += '<table>'
        for line in file:
          line_num += 1
          content  += '<tr><td align="right" valign="top">%d</td><td>%s</td></tr>' \
            % (line_num, saxutils.escape(line))
        content += '</table>'
      except IOError, e:
        content = 'File not found'
      c = RequestContext(request, {'content': mark_safe(content),
                                   'title': file_name,
                                   'pagetitle': file_name})
      t = loader.get_template('file/file_listing.html')
      return HttpResponse(t.render(c))

def members(request):
    t = loader.get_template('main/members.html')
    c = RequestContext(request, {})
    return HttpResponse(t.render(c))

def funding(request):
    t = loader.get_template('main/funding.html')
    c = RequestContext(request, {})
    return HttpResponse(t.render(c))

def research(request):
    t = loader.get_template('main/research.html')
    c = RequestContext(request, {})
    return HttpResponse(t.render(c))

def research_user_interface(request):
    t = loader.get_template('main/research_user_interface.html')
    c = RequestContext(request, {})
    return HttpResponse(t.render(c))

def orthoscope(request, family_accession):
    t = loader.get_template('phyloscope/phyloscope.html')
    c = RequestContext(request, {'json_tree': annotatedTree(family_accession),
                                 'family_accession': family_accession,
                                 'superorthologous_nodes': superorthologousNodes(family_accession),
                                 'hide_subfam_controls': True,
                                 'do_color_ortho': True})
    return HttpResponse(t.render(c))

def orthoscope_get(request):
    t = loader.get_template('phyloscope/phyloscope.html')
    subtreeLeftId = get_parameter_from_GET('subtreeLeftId', request)
    highlightLeftIds = get_parameter_from_GET('highlightLeftIds', request)
    colorLeftIds = get_parameter_from_GET('colorLeftIds', request)
    phog_accession =  get_parameter_from_GET('phog', request)
    family_accession =  get_parameter_from_GET('family', request)
    url_tree_method = get_parameter_from_GET('treeMethod', request)
    if url_tree_method == '':
      url_tree_method = 'ml'
    if phog_accession != "":
      phog_id = int(phog_accession[4:])
      phog = TreeNode.objects.get(id = phog_id)
      if phog == None:
        return render_to_response('404.html')
      else:
        c = RequestContext(request, {
            'json_tree': annotatedTree(phog.tree.family.__unicode__(), 
                            tree_method = url_tree_method),
            'superorthologous_nodes': "",
            'family_accession': phog.tree.alignment.family.__unicode__(),
            'subtree_left_id': phog.left_id,
            'highlight_left_ids': highlightLeftIds,
            'color_left_ids': colorLeftIds,
            'hide_subfam_controls': True,
            'tree_method': url_tree_method,
            'do_color_ortho': True})
        return HttpResponse(t.render(c))
    elif family_accession != "":  
      c = RequestContext(request, {
            'json_tree': annotatedTree(request.GET['family'],
                            tree_method = url_tree_method),
            'superorthologous_nodes': superorthologousNodes(
                      request.GET['family'], tree_method = url_tree_method),
            'family_accession': request.GET['family'],
            'subtree_left_id': subtreeLeftId,
            'highlight_left_ids': highlightLeftIds,
            'color_left_ids': colorLeftIds,
            'hide_subfam_controls': True,
            'do_color_ortho': True,
            'tree_method': url_tree_method,
            'show_family_download': False,
            'show_tree_download': False})
#      add_download_form_to_context(c)
      return HttpResponse(t.render(c))
    else:
      return render_to_response('404.html')

def queryDb(request):
    t = "PhyloFacts::Database Query"
    s = "by accession, taxonomy, biochemistry, localization, process..."

    # Taxon classifications are a bit hacky because we have two that are
    # defined by negation: "invertebrates", which are all metazoans except
    # Chordata, and "prokaryotes", which are all cellular organisms except
    # Eukaryota. To make their designation both meaningful and recognisable,
    # we set it to the additive inverse of the taxon that's not in the group.
    # For instance, Eukaryota is 2759. So, Prokaryota, which consists of
    # Archaea and Bacteria, is -2759 because it's cellular organisms minus
    # eukaryotes.
    taxa = [(5833, 'Plasmodium falciparum'),
            (5782, 'Dictyostelium'),
            (2759, 'Eukaryotes'),
            (7227, 'Fruit fly (Drosophila)'),
            (9606, 'Human (Homo sapiens)'),
            (-7711, 'Invertebrates (metazoans except vertebrates)'),
            (10090, 'Mouse (Mus musculus)'),
            (3193, 'Plants'),
            (-2759, 'Prokaryotes (Archaea and Bacteria)'),
            (7742, 'Vertebrates'),
            (10239, 'Viruses'),
            (4932, 'Yeast (Saccharomyces cerevisiae)')]

    if request.method == 'POST':
	import logging
	logging.debug(request.POST);
	
	# Our goal here is to select sequences. One line for each sequence, fine, but what goes on it?
	# Our sequence_id internal identifier is not appropriate here. UniProt perhaps?

	db = MySQLdb.connect(db=phylofacts_db_name, host=phylofacts_db_host, user=phylofacts_db_user, passwd=phylofacts_db_passwd)
	cur = db.cursor(MySQLdb.cursor.DictCursor)
	
	# Start with set-intersection or union (as the case may be) of families containing target taxa.
	# This is incredibly slow (sometimes >5min) done as a DB join, so we'll get them separately.
	searchTaxa = request.POST.from_taxon + request.POST.containing_sequences_from
	if -7711 in searchTaxa:
	    searchTaxa.remove(-7711)
	    searchTaxa.append(33208)
	if -2759 in searchTaxa:
	    searchTaxa.remove(-2759)
	    searchTaxa.append(2157) # archaea
	    searchTaxa.append(2)    # eubacteria
	sql = """SELECT DISTINCT genbank_uniprot.ncbi_taxid, family.scopid_or_bpgid 
			FROM family, alignment, alignment_protein_sequence, protein_sequence, genbank_uniprot, ncbi_taxonomy
			WHERE family.id = alignment.family_id
			AND alignment.id = alignment_protein_sequence.alignment_id
			AND alignment_protein_sequence.protein_sequence_id = protein_sequence.id
			AND protein_sequence.genbank_uniprot_id = genbank_uniprot.id
			AND genbank_uniprot.ncbi_taxid IN (%s)""" % ','.join(searchTaxa)
	# There's an issue here: some of our categories are quite general, e.g. "eukaryotes", and some have these "excluding"
	# properties, e.g. "invertebrates". We can solve this using left_ids and right_ids in ncbi_taxonomy, e.g.:
	#
	# select distinct genbank_uniprot.ncbi_taxid, family.scopid_or_bpgid from family, alignment, alignment_protein_sequence, 
	# protein_sequence, genbank_uniprot inner join ncbi_taxonomy using (id) inner join ncbi_taxonomy as include on 
	# (ncbi_taxonomy.left_id between include.left_id and include.right_id) where family.id = alignment.family_id and 
	# alignment.id = alignment_protein_sequence.alignment_id and alignment_protein_sequence.protein_sequence_id = protein_sequence.id 
	# and protein_sequence.genbank_uniprot_id = genbank_uniprot.id and include.id = 2157;
	#
	# and, furthermore, using further INNER JOINs to *exclude* stuff ... but how to handle the case where
	# someone wants, oh, I dunno, sequences that also appear in plants and invertebrates?

	return render_to_response('test/query_db.html', {'pagetitle': t, 'subtitle': s, 'taxa': taxa, 'POST': request.POST})

    else:
	form = DbQueryForm()

    return render_to_response('test/query_db.html', {'pagetitle': t, 'subtitle': s, 'form': form, 'taxa': taxa})

def search(request, search_words):
    newPOST = request.POST.copy()
    newPOST.update({'search_words': search_words})
    request.POST = newPOST
    return render_to_response('search/search.html', {'form': SearchForm(request.POST)})
  
def add_download_form_to_context(context):
  family_accession = context['family_accession']
  if family_accession:
    results = Family.objects.filter(scopid_or_bpgid__exact=family_accession)
    if results:
      family_id = results[0].id
      preferred_name = ''
      try:
        preferred_name = context['preferred_name']
      except KeyError:
        pass
      if not preferred_name:
        preferred_name = family_accession
      q = QueryDict('')
      q = q.copy()
      q.update({'family': family_accession})
      q.update({'draft_family': family_accession})
      q.update({'family_id': family_id})
      q.update({'preferred_name': preferred_name})
      q.update({'gathering_method': 'user'})
      q.update({'family_dir': family_dir(family_accession) + 'user'})
      q.update({'workdir': view_tree_work_dir()})
      q.update({'n_subfamilies': '1'}) # Force display of tree download links
      q.update({'n_sequences': '3'})   # (Ditto) Will be superceded below
      q.update({'download_f': '1'})
      view_tree_URL = "view_tree.php?%s" % (q.urlencode())  
      q.update({'view_tree_URL': view_tree_URL})
  
      results = Alignment.objects.filter(family__exact=family_id).filter(subfamily_alignment_id__exact=0)
      if results:
        alignment = results[0]
        q.update({'alignment_id': alignment.id})
        q.update({'seed_sequence_id': alignment.seed_sequence_id})
        q.update({'n_subfamilies': alignment.n_subfamilies}) 
        q.update({'n_sequences': alignment.n_sequences}) 
        q.update({'subfamily': 'Subfamily'})
      context['download_form'] = HiddenDownloadForm(q)
  return

def make_go_summary_line(go_summary):
  return ', '.join(['<a href="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s">%s</a> (<a href="http://www.geneontology.org/GO.evidence.shtml#%s" title="%s">%s</a>)' \
    % (acc, name, evidence.lower(), description, evidence) 
    for (acc, name, evidence, description) in go_summary])
        
class PartnerRow:
  def __init__(self, uniprot, DIP_link_string=None):
    self.uniprot_id = uniprot.__str__().split('_')[0]
    self.accession = uniprot.uniprot_accession
    self.taxon = uniprot.taxon
    self.taxon_url = uniprot.taxon.get_absolute_url()
    self.description = uniprot.description
    self.go_experimental_summaries = {}
    self.go_nonexperimental_summaries = {}
    self.ppi_link = uniprot.get_ppi_link()
    self.DIP_link = DIP_link_string
    for go_component in ['biological_process', 'molecular_function',
                          'cellular_component']:
      long_summary = []
      for experimental in [True, False]:
        summary = uniprot.get_go_summary(go_component, experimental=experimental)
        if summary:
          long_summary = long_summary + summary
      short_summary = []
      current_length = 0
      for (acc, name, evidence, description) in long_summary:
        next_str = ", %s (%s)" % (name, evidence)
        if len(next_str) + current_length <= 200:
          short_summary = short_summary + \
                            [(acc, name, evidence, description)]
          current_length += len(next_str)
        else:
          break
      attr_name = "go_short_summary_%s" % go_component
      setattr(self, attr_name, mark_safe(make_go_summary_line(short_summary)))
      if len(short_summary) != len(long_summary):
        attr_name = "go_long_summary_%s" % go_component
        setattr(self, attr_name, mark_safe(make_go_summary_line(long_summary)))
        
def orthologs_about(request):
  context = {}
  return render_to_response('phog/about_phog.html', context)

def orthologs_assessment(request):
  context = {}
  return render_to_response('phog/phog_supplement_2.html', context)

def orthologs_downloads(request):
    if 'taxon1' not in request.GET or 'taxon2' not in request.GET:
      return render_to_response('phog/downloads.html',
                              {'form': PairwiseGenomicOrthologyForm()})
    else:
      form = PairwiseGenomicOrthologyForm({'taxon1': request.GET['taxon1'],
                                          'taxon2': request.GET['taxon2']})
      if not form.is_valid():
        return render_to_response('phog/downloads.html', {'form': form})
      else:
        taxon_id1 = int(form.cleaned_data['taxon1'])
        taxon_id2 = int(form.cleaned_data['taxon2'])
        if taxon_id1 > taxon_id2:
          # Switch taxon_id1 and taxon_id2
          taxon_id2 += taxon_id1
          taxon_id1 = taxon_id2 - taxon_id1
          taxon_id2 -= taxon_id1
        taxon_object1 = UniProtTaxonomy.objects.get(id = taxon_id1)
        taxon_object2 = UniProtTaxonomy.objects.get(id = taxon_id2)
        taxon_name1 = '_'.join(taxon_object1.scientific_name.split()[0:2])
        taxon_name2 = '_'.join(taxon_object2.scientific_name.split()[0:2])
        filepattern ='/clusterfs/vasudha/bpg/OrthologsForQuest/orthoXML/%s-%s_*_orthologs.xml' % (taxon_name1, taxon_name2)
        files = glob.glob(filepattern)
        if len(files) == 0:
          return render_to_response('phog/downloads.html', 
                      {'form': form, 'missing_file_error': 
          "The pairwise genomic ortholog file between %s and %s is not present"
                      % (taxon_name1, taxon_name2)})
        file = files[0]
        filename = os.path.split(file)[1]
        f = open(file)
        response = HttpResponse(mimetype="text/xml")
        response['Content-Disposition'] = 'attachment; filename=%s' % filename
        response.write(f.read())
        f.close()
        return response

def orthologs_faq(request):
  context = {}
  return render_to_response('phog/phog_faq.html', context)

def orthologs_quickstart(request):
    context = { }
    return render_to_response('phog/phog_quick_start_guide.html', context)

def orthologs_supplement(request, version='v1'):
  context = {}
  return render_to_response('phog/phog_supplement_1.html', context)

orthologs_tutorial_templates = [
    'tutorial-01.html',
    'tutorial-02.html',
    'tutorial-03.html',
    'tutorial-04.html',
    'tutorial-05.html',
    'tutorial-06.html',
    'tutorial-07.html',
    'tutorial-08.html',
    'tutorial-09.html',
    'tutorial-10.html',
    'tutorial-11.html',
    'tutorial-12.html',
    'tutorial-13.html',
    'tutorial-14.html',
    'tutorial-15.html',
    'tutorial-16.html',
    'tutorial-17.html',
    'tutorial-18.html',
    'tutorial-19.html',
    'tutorial-20.html'
    ]
 
def orthologs_tutorial(request):
    page_num = 1
    if 'page' in request.GET:
        try:
            page_num = int(request.GET['page'])
        except ValueError:
            return render_to_response('404.html')
        
    context = { 'page' : page_num }
    context['padded_page_num'] = '%02d' %page_num
    context['page_numbers'] = range(1, len(orthologs_tutorial_templates) + 1)
    
    if page_num > 1 :
        context['previous'] = page_num - 1
        context['padded_previous'] = '%02d' %(page_num - 1)
    if page_num < len(orthologs_tutorial_templates) :
        context['next'] = page_num + 1
        context['padded_next'] = '%02d' %(page_num + 1)
    
    return render_to_response('phog/tutorial/%s'
                              % orthologs_tutorial_templates[page_num  - 1],
                              context)

def crispr(request):
  context = {'page_title': 'CRISPR related restriction enzyme families'}
  context_for_datagrid = {}
  context_for_datagrid['extra_query'] = ''
  families = Family.objects.filter(id__in = 
    [81908, 81909, 81916, 81905, 81898, 81907, 81900, 81899, 81910, 81904,
    81902, 81906, 81911, 81920, 81914, 81915, 81927, 81918, 81917, 81919,
    81922, 81928, 81923, 81926, 81930, 81929, 81897, 81901, 81903, 81912,
    81913, 81921, 81925, 81924, 81931, 81953, 81955, 81954, 81952, 81956,
    81957, 81958, 81959, 81934, 81932, 81935, 81936, 81937, 81933, 81938,
    81940, 81941, 81943, 81946, 81945, 81944, 81942, 81939, 81947, 81948])
  return FamilyDataGrid(request, families, 
        extra_context=context_for_datagrid).render_to_response(
                                    'report/crispr.html', context)

def phog(request, phog_accession=None):
  if 'phog_accession' in request.GET:
    putative_phog_accession = request.GET['phog_accession'].rstrip()
    if len(putative_phog_accession) == 6:
      m = uniprot_accession_re1.match(putative_phog_accession.upper())
      if not m:
        m = uniprot_accession_re2.match(putative_phog_accession.upper())
      if m and len(m.group()) == len(putative_phog_accession):
        redirect_url = "/phog/orthologs/%s" % putative_phog_accession
        return HttpResponseRedirect(redirect_url)
    m = uniprot_identifier_re1.match(putative_phog_accession.upper())
    if not m:
      m = uniprot_identifier_re2.match(putative_phog_accession.upper())
    if not m:
      m = uniprot_identifier_re3.match(putative_phog_accession.upper())
    if not m:
      m = gi_re.match(putative_phog_accession)
    if m and len(m.group()) == len(putative_phog_accession):
      redirect_url = "/phog/orthologs/%s" % putative_phog_accession
      return HttpResponseRedirect(redirect_url)
    m = scopid_re.match(putative_phog_accession.lower())
    if not m:
      m = bpgid_re.match(putative_phog_accession.lower())
    if m and len(m.group()) == len(putative_phog_accession):
      redirect_url = "/family/%s" % putative_phog_accession
      return HttpResponseRedirect(redirect_url)
  context = {'page_title': "PHOG: PhyloFacts Orthology Group"}
  context_for_datagrid = {}
  if phog_accession != None or \
      ('phog_accession' in request.GET and \
        request.GET['phog_accession'] != request.GET['phog_accession'].upper()):
    if not phog_accession:
      phog_accession = request.GET['phog_accession'].rstrip()
    phog_accession = phog_accession.upper().rstrip()
    newGET = request.GET.copy()
    newGET.update({'phog_accession': phog_accession})
    request.GET = newGET
    context['phog_accession'] = phog_accession
    context_for_datagrid['extra_query'] = "phog_accession=%s&" % phog_accession
    form = PHOGForm(request.GET)
    form.phog_accession = phog_accession
  elif 'phog_accession' in request.GET:
    phog_accession = request.GET['phog_accession'].rstrip()
    context['phog_accession'] = request.GET['phog_accession']
    context_for_datagrid['extra_query'] \
      = "phog_accession=%s&" % request.GET['phog_accession']
    form = PHOGForm(request.GET)
    form.phog_accession = request.GET['phog_accession']
  else:
    form = PHOGForm()
  context['form'] = form
  if form.is_bound and form.is_valid() and phog_accession != None and \
      phog_accession != u'':
    ortholog_type, threshold \
      = get_ortholog_type_threshold_from_phog_accession(phog_accession)
    (phog, orthologs, error) \
      = getPHOGQuerySet(form.cleaned_data['phog_accession'], 
                        ortholog_type = ortholog_type,
                        threshold = threshold)
    if (error != None):
        context['db_error'] = error
    else:
        phog_row = PHOGRow(phog, ortholog_type, threshold)
        context['phog'] = phog_row
        return OrthologDataGrid(request, orthologs, 
              extra_context=context_for_datagrid,
              includePHOGColumns=False).render_to_response(
                                          'phog/phog.html', context)
  return render_to_response('phog/phog.html', context)

def interactors(request, sequence_id=None):
    if 'sequence_id' in request.GET:
      putative_sequence_id = request.GET['sequence_id']
      m = phog_re.match(putative_sequence_id.upper())
      if m and len(m.group()) == len(putative_sequence_id):
        redirect_url = "/phog/%s" % putative_sequence_id
        return HttpResponseRedirect(redirect_url)
      m = scopid_re.match(putative_sequence_id.lower())
      if not m:
        m = bpgid_re.match(putative_sequence_id.lower())
      if m and len(m.group()) == len(putative_sequence_id):
        redirect_url = "/family/%s" % putative_sequence_id
        return HttpResponseRedirect(redirect_url)
    context = {'page_title': "Protein-Protein Interactions"}
    context_for_datagrid = {}
    if sequence_id == None and \
          ('sequence_id' not in request.GET or \
            request.GET['sequence_id'] == u'') and \
        'sequence_fasta' in request.GET and \
        request.GET['sequence_fasta'] != u'':
      m = fasta_re.match(request.GET['sequence_fasta'])
      if m:
        (sequence_id, record) = get_uniprot_id_from_fasta(request.GET['sequence_fasta'])
        if sequence_id == None:
          context['fasta_error'] \
            = 'No known interactions with your sequence.'
      else:
        context['form_error'] = 'Invalid FASTA sequence'
    if sequence_id != None or \
        ('sequence_id' in request.GET and \
          request.GET['sequence_id'] != request.GET['sequence_id'].upper()):
      if not sequence_id:
        sequence_id = request.GET['sequence_id']
      sequence_id = sequence_id.upper()
      newGET = request.GET.copy()
      newGET.update({'sequence_id': sequence_id})
      request.GET = newGET
      context['sequence_id'] = sequence_id
      context_for_datagrid['extra_query'] = "sequence_id=%s&" % sequence_id
      form = OrthologForm(request.GET)
      form.sequence_id = sequence_id
    elif 'sequence_id' in request.GET:
      sequence_id = request.GET['sequence_id']
      context['sequence_id'] = request.GET['sequence_id']
      context_for_datagrid['extra_query'] \
        = "sequence_id=%s&" % request.GET['sequence_id']
      form = OrthologForm(request.GET)
      form.sequence_id = request.GET['sequence_id']
    else:
      form = OrthologForm()
    context['form'] = form
    if form.is_bound and form.is_valid() and sequence_id != None and \
        sequence_id != u'':
      if form.non_field_errors():
        context['form_error'] = form.non_field_errors().as_text()
      else:
        (uniprot, uniprots, error) = getInteractorQuerySet(sequence_id)
        if error:
          context['db_error'] = error
        else:
          context['query_description'] = uniprot.de
          context['query_species'] = uniprot.taxon.scientific_name
          return PartnerDataGrid(request, uniprots, query_uniprot=uniprot,
                  extra_context=context_for_datagrid).render_to_response(
                    'phog/interactors.html', context)
    return render_to_response('phog/interactors.html', context)


def orthologs(request, sequence_id=None):
    with_family_type = request.GET.get('with_family_type', False)
    for_iframe = request.GET.get('for_iframe', False)
    if 'sequence_id' in request.GET:
      putative_sequence_id = request.GET['sequence_id'].strip()
      m = phog_re.match(putative_sequence_id.upper())
      if m and len(m.group()) == len(putative_sequence_id):
        redirect_url = "/phog/%s" % putative_sequence_id
        return HttpResponseRedirect(redirect_url)
      m = scopid_re.match(putative_sequence_id.lower())
      if not m:
        m = bpgid_re.match(putative_sequence_id.lower())
      if m and len(m.group()) == len(putative_sequence_id):
        redirect_url = "/family/%s" % putative_sequence_id
        return HttpResponseRedirect(redirect_url)
    context = {'page_title': "PHOGs: PhyloFacts Orthology Groups"}
    context_for_datagrid = {}
    if sequence_id == None and \
          ('sequence_id' not in request.GET or \
            request.GET['sequence_id'] == u'') and \
        'sequence_fasta' in request.GET and \
        request.GET['sequence_fasta'] != u'':
      m = fasta_re.match(request.GET['sequence_fasta'])
      if m:
        (sequence_id, record) = get_uniprot_id_from_fasta(request.GET['sequence_fasta'])
        sequence_id = None
        if not record is None:
          return approximate_matches(request, record, context)
        else:
          context['fasta_error'] = 'Invalid FASTA sequence'
      else:
        context['fasta_error'] = 'Invalid FASTA sequence'
    if sequence_id != None or \
        ('sequence_id' in request.GET and \
          request.GET['sequence_id'] != request.GET['sequence_id'].upper()):
      if not sequence_id:
        sequence_id = request.GET['sequence_id'].rstrip()
      sequence_id = sequence_id.upper().rstrip()
      newGET = request.GET.copy()
      newGET.update({'sequence_id': sequence_id})
      request.GET = newGET
      context['sequence_id'] = sequence_id
      context_for_datagrid['extra_query'] = "sequence_id=%s&" % sequence_id
      form = OrthologForm(request.GET)
      form.sequence_id = sequence_id
    elif 'sequence_id' in request.GET:
      sequence_id = request.GET['sequence_id'].rstrip()
      context['sequence_id'] = request.GET['sequence_id']
      context_for_datagrid['extra_query'] \
        = "sequence_id=%s" % request.GET['sequence_id']
      form = OrthologForm(request.GET)
      form.sequence_id = request.GET['sequence_id']
    else:
      form = OrthologForm()
    context['form'] = form
    if form.is_bound and form.is_valid() and sequence_id != None and \
        sequence_id != u'':
      if form.non_field_errors():
        context['form_error'] = form.non_field_errors().as_text()
      else:
        # For some reason, ortholog_type comes in as Unicode
        try:
          ortholog_type = int(form.cleaned_data['ortholog_type']) 
        except ValueError:
          ortholog_type = OrthologTypes.SuperOrtholog
        threshold = form.cleaned_data['threshold']
        context_for_datagrid['extra_query'] \
          = context_for_datagrid['extra_query'] \
          + "&ortholog_type=%d&threshold=%g" % (ortholog_type, threshold)
        if with_family_type:
            context_for_datagrid['extra_query'] += '&with_family_type=1'
        if for_iframe:
            context_for_datagrid['extra_query'] += '&for_iframe=1'
        (ownUniProt, phogs, best_phogs, 
            orthologs, phog_of_ortholog, error) \
          = getOrthologQuerySet(form.cleaned_data['sequence_id'],
                                ortholog_type, threshold)
        if (error != None):
            context['db_error'] = error
        else:
            have_nontrivial_phogs = False
            for phog in phogs.all():
              if phog.get_num_nonredundant_sequences(ortholog_type,
                                                      threshold) >= 2:
                have_nontrivial_phogs = True
            if not have_nontrivial_phogs:
              context['db_error'] = not_in_phog(sequence_id)
              return render_to_response('phog/orthologs.html', context)
            phog_rows = [PHOGRow(phog, ortholog_type, threshold) 
                          for phog in phogs]
            context['phogs'] = phog_rows
            
            #context['phogs'] = []
            context['query_description'] = ownUniProt.de
            context['query_species'] = ownUniProt.taxon.scientific_name

            # Determine if the ppi link should be shown and set appropriate
            # context variable if it should.
            context['show_ppi_link'] = False
            # FIXME: 2010/06/13 RSD
            # Commenting this out since we currently don't have the PPI data in
            # the schema or models in the database, let alone populated.
            """
            for ortholog in orthologs:
              if bool(ortholog.sequence_header.uniprot.get_interacting_partners()):
                context['show_ppi_link'] = True
                # context['ortholog_type'] = ortholog_type
                context['netscope_url'] = \
                 ortholog.sequence_header.uniprot.get_netscope_url(ortholog_type, threshold)
                break
            """

            ## @@@ this is not sufficient as the sequence header identifier may differ
            ##     from the sequence id, resulting in no url being stored.
            # Put url for the query into a context variable.
            # for ortholog in orthologs:
            #   if sequence_id == ortholog.sequence_header.identifier():
            #     context['query_url'] = ortholog.sequence_header.get_absolute_url()

            n_sequences_of_family = {}
            for phog in phogs:
              n_sequences = phog.tree.family.canonical_root_node(
                  ).get_num_included_leaves(OrthologTypes.PHOG_T_Custom, 
                                            threshold = 10000.0)
              n_sequences_of_family[phog.tree.family] = n_sequences
            return OrthologDataGrid(request, orthologs, 
                  with_family_type = with_family_type,
                  uniprot = ownUniProt,
                  ortholog_type = ortholog_type,
                  threshold = threshold,
                  phog_of_ortholog = phog_of_ortholog,
                  n_sequences_of_family = n_sequences_of_family,
                  extra_context=context_for_datagrid).render_to_response(
                                          for_iframe and 'phog/orthologs_for_iframe.html' or 'phog/orthologs.html',
                                          context)
    return render_to_response('phog/orthologs.html', context)

def approximate_matches(request, record, context):
    context = {'page_title': "Approximately Matching Sequences"}
    context_for_datagrid = {}
    form = OrthologForm(request.GET)
    context['form'] = form
    extra_context = {}
    (orthologs, alignment_dict, error) \
      = get_approximately_matching_sequence_query_set(record, context, extra_context)
    if (error != None):
      context['db_error'] = error
      return render_to_response('phog/orthologs.html', context)
    else:
      context['approximate_matches'] = True
      return ApproximateMatchesDataGrid(request, orthologs, alignment_dict, extra_context).render_to_response( 
                                        'phog/orthologs.html', context)

def YN_of_bool(flag):
  if flag:
    return 'Y'
  else:
    return 'N'


def csv_ortholog_table_response(ortholog_set, error, filename, 
                                includePHOG = True, ortholog_type =
                                OrthologTypes.SuperOrtholog,
                                threshold = 0.0,
                                phog_of_ortholog = {}):
  response = HttpResponse(mimetype='text/csv')
  response['Content-Disposition'] = 'attachment; filename=%s' % filename
  if (error == None):
    writer = csv.writer(response)
    row = ['Gene','Species','Description','In_SwissProt?',
                'Has_Experimental_Evidence?', 'EC','KEGG','PFAM',
                'Has_Literature?',
                'Family','Global_Homology?']
    if includePHOG:
      row = row + ['PHOG']
    writer.writerow(row)
    for ortholog in ortholog_set:
      sequence_header = ortholog.sequence_header
      pfams = [pfam_tuple[0].name for pfam_tuple in
                    ortholog.tree.family.get_pfams()]
      pfam_str = '"%s"' % ','.join(pfams)
      if sequence_header.uniprot:
        row=[sequence_header.identifier(),
              sequence_header.uniprot.taxon.__str__().replace(' ', '_'),
              sequence_header.description(),
              YN_of_bool(sequence_header.uniprot.in_swissprot_f == 1),
              YN_of_bool(sequence_header.uniprot.has_experimental_evidence()),
              sequence_header.uniprot.get_ec(),
              ','.join(sequence_header.uniprot.get_kegg_map_ids()),
              pfam_str,
              YN_of_bool(sequence_header.uniprot.has_literature()),
              ortholog.tree.family.__str__(),
              YN_of_bool(ortholog.tree.family.is_global_homology())]
      else:
        row=[sequence_header.identifier(),
              '',
              sequence_header.description(),
              '',
              '',
              '',
              '',
              pfam_str,
              '',
              ortholog.tree.family.__str__(),
              YN_of_bool(ortholog.tree.family.is_global_homology())]
      if includePHOG and ortholog in phog_of_ortholog:
        phog = phog_of_ortholog[ortholog]
        row = row + [phog.get_accession(ortholog_type = ortholog_type,
                                        threshold = threshold)]
      writer.writerow(row)
      
  return response
  
def phog_as_csv(request, phog_accession):
  ortholog_type, threshold \
    = get_ortholog_type_threshold_from_phog_accession(phog_accession)
  (phog, ortholog_set, error) = getPHOGQuerySet(phog_accession, 
                                          ortholog_type = ortholog_type, 
                                          threshold = threshold)
  filename = "%s.csv" % phog_accession
  return csv_ortholog_table_response(ortholog_set, error, filename,
                                      includePHOG = False)

def hyperscope(request, phog_accession, taxon_id = None):
  phog_id = int(phog_accession[4:13])
  phog = TreeNode.objects.get(id = phog_id)
  if phog == None:
    return render_to_response('404.html')
  else:
    t = loader.get_template('ppi/hyperscope.html')
    ortholog_type, threshold \
        = get_ortholog_type_threshold_from_phog_accession(phog_accession)
    pig_res = run_interactome_viewer(phog, ortholog_type, threshold, taxon_id)
    if pig_res:
      c = RequestContext(request,
                         {'json_hyper': hyperNeighborhood(phog, ortholog_type,
                         threshold, taxon_id),
                          'node_coords': pig_res[0],
                          'filename': pig_res[1],
        'taxon': taxon_id})
    else:
      return render_to_response('404.html')
  return HttpResponse(t.render(c))

# # Much of this is copied from orthologs.
def netscope(request, sequence_id=None):
    if 'sequence_id' in request.GET:
      putative_sequence_id = request.GET['sequence_id'].rstrip()
      m = phog_re.match(putative_sequence_id.upper())
      if m and len(m.group()) == len(putative_sequence_id):
        redirect_url = "/net/%s" % putative_sequence_id
        return HttpResponseRedirect(redirect_url)
      #m = scopid_re.match(putative_sequence_id.lower())
      #if not m:
      #  m = bpgid_re.match(putative_sequence_id.lower())
      #if m and len(m.group()) == len(putative_sequence_id):
      #  redirect_url = "/family/%s" % putative_sequence_id
      #  return HttpResponseRedirect(redirect_url)
    context = {'page_title': "PHOG Interactome Viewer"}
    context_for_datagrid = {}
    #if sequence_id == None and \
    #      ('sequence_id' not in request.GET or \
    #        request.GET['sequence_id'] == u'') and \
    #    'sequence_fasta' in request.GET and \
    #    request.GET['sequence_fasta'] != u'':
    #  m = fasta_re.match(request.GET['sequence_fasta'])
    #  if m:
    #    (sequence_id, record) = get_uniprot_id_from_fasta(request.GET['sequence_fasta'])
    #    sequence_id = None
    #    if not record is None:
    #      return approximate_matches(request, record, context)
    #    else:
    #      context['fasta_error'] = 'Invalid FASTA sequence'
    #  else:
    #    context['fasta_error'] = 'Invalid FASTA sequence'
    if sequence_id != None or \
        ('sequence_id' in request.GET and \
          request.GET['sequence_id'] != request.GET['sequence_id'].upper()):
      if not sequence_id:
        sequence_id = request.GET['sequence_id'].rstrip()
      sequence_id = sequence_id.upper().rstrip()
      newGET = request.GET.copy()
      newGET.update({'sequence_id': sequence_id})
      request.GET = newGET
      context['sequence_id'] = sequence_id
      context_for_datagrid['extra_query'] = "sequence_id=%s&" % sequence_id
      form = OrthologForm(request.GET)
      form.sequence_id = sequence_id
    elif 'sequence_id' in request.GET:
      sequence_id = request.GET['sequence_id'].rstrip()
      context['sequence_id'] = request.GET['sequence_id']
      context_for_datagrid['extra_query'] \
        = "sequence_id=%s" % request.GET['sequence_id']
      form = OrthologForm(request.GET)
      form.sequence_id = request.GET['sequence_id']
    else: 
      form = OrthologForm()
    context['form'] = form
    if form.is_bound and form.is_valid() and sequence_id != None and \
        sequence_id != u'':
      if form.non_field_errors():
        context['form_error'] = form.non_field_errors().as_text()
      else:
        # For some reason, ortholog_type comes in as Unicode
        try:
          ortholog_type = int(form.cleaned_data['ortholog_type']) 
        except ValueError:
          ortholog_type = OrthologTypes.SuperOrtholog
        threshold = form.cleaned_data['threshold']
        # threshold = preset_thresholds[ortholog_type]
        if (OrthologTypes.SuperOrtholog == ortholog_type):
          context['threshold_type']  = "Super orthology"
        elif (OrthologTypes.PHOG_T_Tight == ortholog_type):
          context['threshold_type']  = "Close"
        elif (OrthologTypes.PHOG_T_Medium == ortholog_type):
          context['threshold_type']  = "Moderate"
        elif (OrthologTypes.PHOG_T_Loose == ortholog_type):
          context['threshold_type']  = "Distant"
        elif (OrthologTypes.PHOG_T_Custom == ortholog_type):
          context['threshold_type']  = "Custom"
        else:
          context['threshold_type']  = "Unknown"
        context_for_datagrid['extra_query'] \
          = context_for_datagrid['extra_query'] \
          + "&ortholog_type=%d&threshold=%g" % (ortholog_type, threshold)
        (ownUniProt, phogs, best_phogs, 
            orthologs, phog_of_ortholog, error) \
          = getOrthologQuerySet(form.cleaned_data['sequence_id'],
                                ortholog_type, threshold)
        if (error != None):
            context['db_error'] = error
        else:
            have_nontrivial_phogs = False
            for phog in phogs.all():
              if phog.get_num_nonredundant_sequences(ortholog_type,
                                                      threshold) >= 2:
                have_nontrivial_phogs = True
            if not have_nontrivial_phogs:
              context['db_error'] = not_in_phog(sequence_id)
              return render_to_response('ppi/netscope.html', context)
            context['query_description'] = ownUniProt.description
            context['query_species'] = ownUniProt.taxon.scientific_name
            
            taxon = ownUniProt.taxon
            taxon_id = taxon.id

            (edge_list, node_coords, filename, observed_partners, predicted_partners) = \
                        run_netscope(ownUniProt, taxon, ortholog_type, threshold)
            context['edge_list']         = edge_list
            context['node_coords']       = node_coords
            context['filename']          = filename
            context['taxon_scientific']  = taxon.scientific_name
            context['taxon_common']      = taxon.common_name
            context['query_sequence_id'] = ownUniProt.id
            context['threshold']         = threshold
            context['ortholog_type']     = ortholog_type

            observed_partner_rows = []
            for partner in observed_partners:
              observed_partner_rows.append(PartnerRow(partner[0], partner[1]))
            context['observed_partners'] = observed_partner_rows

            predicted_partner_rows = []
            for partner in predicted_partners:
              predicted_partner_rows.append(PartnerRow(partner))
            context['predicted_partners'] = predicted_partner_rows

            return render_to_response('ppi/netscope.html', context)
    return render_to_response('ppi/netscope.html', context)


    
def phog2phog(request, phog_accession1, phog_accession2):
  phog_id1 = int(phog_accession1[4:13])
  phog1 = TreeNode.objects.get(id = phog_id1)
  phog_id2 = int(phog_accession2[4:13])
  phog2 = TreeNode.objects.get(id = phog_id2)
  ortholog_type, threshold \
    = get_ortholog_type_threshold_from_phog_accession(phog_accession1)
  ortholog_type2, threshold2 \
    = get_ortholog_type_threshold_from_phog_accession(phog_accession2)
  if phog1 == None or phog2 == None or ortholog_type != ortholog_type2 \
      or threshold != threshold2:
    return render_to_response('404.html')
  else:
    phog_row1 = PHOGRow(phog1, ortholog_type, threshold)
    phog_row2 = PHOGRow(phog2, ortholog_type, threshold)
    context = { 'phog_row1': phog_row1, 'phog_row2': phog_row2 }
    edges, uniprot1_of_edge, uniprot2_of_edge \
        = getPHOG2PHOGQuerySet(phog1, phog2, ortholog_type, threshold)
    return PHOG2PHOGDataGrid(request, edges, uniprot1_of_edge,
                              uniprot2_of_edge).render_to_response(
                                  'phog/phog2phog.html', context)

def phog_tree(request, phog_accession, leftid1=None, leftid2=None):
  tree_id = int(phog_accession[4:11])
  phog_left_id = int(phog_accession[12:])
  phog = TreeNode.objects.get(tree__id = tree_id, left_id = phog_left_id)
  if phog == None:
    return render_to_response('404.html')
  else:
    t = loader.get_template('phyloscope/phyloscope.html')
    highlightLeftIds = ""
    if leftid1 or leftid2:
      if not leftid2:
        highlightLeftIds = leftid1
      elif not leftid1:
        highlightLeftIds = leftid2
      else:
        highlightLeftIds = "%s %s" % (leftid1, leftid2)
    colorLeftIds = ''
    c = RequestContext(request, {
                    'json_tree': annotatedTree(phog.tree.family.__unicode__(),
                    tree_method=phog.tree.method),
                    'superorthologous_nodes': "",
                    'family_accession': phog.tree.family.__unicode__(),
                    'phog_accession': phog_accession,
                    'subtree_left_id': phog.left_id,
                    'highlight_left_ids': highlightLeftIds,
                    'color_left_ids': colorLeftIds,
                    'hide_subfam_controls': True,
                    'hide_ortho_controls': True,
                    'show_full_tree_download': False})
#    add_download_form_to_context(c)
    return HttpResponse(t.render(c))

def phog_tree_json(request, phog_accession):
  tree_id = int(phog_accession[4:11])
  phog_left_id = int(phog_accession[12:])
  phog = TreeNode.objects.get(tree__id = tree_id, left_id = phog_left_id)
  if phog == None:
    return render_to_response('404.html')
  response = HttpResponse(mimetype='application/json')
  response['Content-Disposition'] = 'inline; filename=%s.json' % phog_accession
  response.write(annotatedTree(phog.tree.family.__unicode__(), 
                  tree_method=phog.tree.method))
  return response

def interactors_as_csv(request, sequence_id):
  (interactor_set, error) = getInteractorSet(sequence_id)
  filename = "%s-interactors.csv" % sequence_id
  response = HttpResponse(mimetype='text/csv')
  response['Content-Disposition'] = 'attachment; filename=%s' % filename
  if (error == None):
    writer = csv.writer(response)
    writer.writerow(['UniProt_ID','Description'])
    for interactor in interactor_set:
      writer.writerow([interactor.uniprot_id, interactor.description()])
  return response

def orthologs_as_csv(request, sequence_id):
  ortholog_type = OrthologTypes.SuperOrtholog
  threshold = 0.0
  if 'ortholog_type' in request.GET:
    ortholog_type_str = request.GET['ortholog_type']
    m = gi_re.match(ortholog_type_str)
    if m and len(m.group()) == len(ortholog_type_str):
      try:
        ortholog_type = int(ortholog_type_str)
        if ortholog_type < OrthologTypes.SuperOrtholog \
            or ortholog_type > OrthologTypes.PHOG_T_Custom:
          ortholog_type = OrthologTypes.SuperOrtholog
      except ValueError:
        ortholog_type = OrthologTypes.SuperOrtholog
    if ortholog_type == OrthologTypes.PHOG_T_Custom and \
        'threshold' in request.GET:
      threshold_str = request.GET['threshold']
      m = float_re.search(threshold_str)
      if m and len(m.group()) == len(threshold_str):
        try:
          threshold = float(threshold_str)
          if threshold < 0.0:
            threshold = 0.0
        except ValueError:
          threshold = 0.0
  (ownUniProt, phogs, best_phogs, 
   ortholog_set, phog_of_ortholog, error) \
    = getOrthologQuerySet(sequence_id, ortholog_type = ortholog_type,
                          threshold = threshold)
  if ortholog_type == OrthologTypes.SuperOrtholog:
    suffix = ''
  elif ortholog_type == OrthologTypes.PHOG_T_Tight:
    suffix = 'TC'
  elif ortholog_type == OrthologTypes.PHOG_T_Medium:
    suffix = 'TM'
  elif ortholog_type == OrthologTypes.PHOG_T_Loose:
    suffix = 'TD'
  elif ortholog_type == OrthologTypes.PHOG_T_Custom:
    suffix = 'T%g' % threshold
  filename = "%s-orthologs_%g.csv" % (sequence_id, threshold)
  return csv_ortholog_table_response(ortholog_set, error, filename, 
                                      includePHOG = True, 
                                      ortholog_type = ortholog_type,
                                      threshold = threshold,
                                      phog_of_ortholog = phog_of_ortholog)

def home(request):
  redirect_url = "/phylofacts/"
  return HttpResponseRedirect(redirect_url)

def family(request, family_accession):
  redirect_url \
    = "/phylofacts/family/family_info.php?family=%s" \
        % family_accession
  return HttpResponseRedirect(redirect_url)

def intrepid(request):
  redirect_url = "/intrepid/"
  return HttpResponseRedirect(redirect_url)

