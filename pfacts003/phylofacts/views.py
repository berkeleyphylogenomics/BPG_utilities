import os, re, cPickle, string, logging, glob
import sys
import time
if sys.version_info < (2,5):
    import simplejson
else:
    import json as simplejson
import pfacts003.utils.json as json
from array import array
from Bio import AlignIO, SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.pairwise2 import dictionary_match
from django.shortcuts import render_to_response, redirect
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.db import connection
import hashlib, random
from pfacts003.phylofacts.models import Family, OrthologTypes, UniProt, \
    TreeNode, SequenceHeader, TreeNodeConsensus, HMM, PDB, PHOGRow, FamilyCache, \
    UniProtTaxonomy, Tree, Pfam, MetacycNode, FatcatJob, FatcatJobFamily, FxnSitePredictionJob, MetacycReaction, OhanaQueueStatus, OhanaJobs
from pfacts003.phylofacts.models import \
    get_ortholog_type_threshold_from_phog_accession
from pfacts003.phylofacts.datagrids import FamilyDataGrid
from pfacts003.phog.orthologs import getOrthologQuerySet, getPHOGQuerySet
from pfacts003.phog.datagrids import OrthologDataGrid
from pfacts003.phog.forms import OrthologForm
from pfacts003.phog.views import PHOGRow
from pfacts003.utils.json_data_objects import get_alignment_rows, \
  annotatedSummaryTree
from pfacts003.utils.get_identifier_from_fasta import \
                 get_uniprot_id_from_fasta, \
                 get_approximately_matching_family_query_set, \
                 get_seq_record_from_defline_and_sequence, \
                 get_single_uniprot_object_from_seq_record, \
                 get_family_data_for_sequence
from pfacts003.utils.intrepid import parse_intrepid_output
from pfacts003.utils.id_patterns import fasta_re, is_uniprot_identifier_format, is_uniprot_accession_format, bpgid_re, phog_re
from pfacts003.utils.links import make_pfam_links, alignment_hmm_pairwise_link
from django.utils.safestring import mark_safe
from django.template import Context, RequestContext, loader
from pfacts003.utils.files import PHOG_BLAST_DIR
from urllib import unquote_plus
from pfacts003.utils.consts import max_phyloscope_tree_size
from pfacts003.phylofacts.forms import SequenceSearchForm, FamilyAccessionForm, UniProtForm, NewSequenceSearchForm
from pfacts003.phylofacts.models import get_dir_of_family_accession
from pfacts003.utils import messages
from pfacts003.utils.extract_sequence_from_fasta import extract_sequence_from_fasta
from pfacts003.utils.annotation import summarize
from pfacts003.utils.bpg_queue import ohana_queue_status
from textwrap import wrap

logging.basicConfig(level = logging.ERROR)

trivial_translation = string.maketrans('', '')
dotlowercase = '.' + string.lowercase
dotdash = '.-'
contiguous_lowercase = re.compile('[a-z]+')

blosum62_of_residues = dictionary_match(SubsMat.SeqMat(MatrixInfo.blosum62))

wrapwidth = 80

def get_parameter_from_GET(parameter, request):
    if parameter in request.GET:
        return request.GET[parameter]
    else:
        return ""

# TODO: this should be removed. The code is needlessly dense. Also,
# rendering a 404 is only one possible response to not finding a UniProt
# record in our db; we need to allow a more graceful way of handling this as well.
# This def should be replaced by
# phylofacts.models.UniProt.find_by_accession_or_identifier() classmethod.
# comment by ST,
def get_up_by_acc_or_ident_or_404(acc_or_ident):
    '''what a verbose name'''
    try:
        return UniProt.objects.get(**{
            len(acc_or_ident) > 6 and 'uniprot_identifier' \
            or 'accession': acc_or_ident
        })
    except UniProt.DoesNotExist:
        raise Http404

def search(request):
    """Solr Search.
    This function now just redirects to a separate search page for the solr search.
    """
    if not ( 'query' in request.GET ):
        raise Http404

    query = request.GET['query']

    return render_to_response('phylofacts/search.html', { 'query': query })


def accession_search(request):
    """Accession Search

    A search field is placed in the upper-right hand corner of all pages. There are also fields to 'Jump to a PhyloFacts family' and to enter an accession for searching.

    action of that search is directed to this view via a GET method.

    Currently, this simply redirects the search value to the UniProt page
    /phylofacts/sequence/UniProt.
    """

    query = request.GET.get('query', None)
    search_key = request.GET.get('search_key', None)

    if search_key == 'bpg' or query.startswith('bpg'):
        if not re.match('^bpg\d{7}$', query):
            response = HttpResponse()
            response.write("Bad family format. Should be: bpgXXXXXXX <br />\n")
            return response
        else:
            return HttpResponseRedirect('/phylofacts/family/%s/'% query)

    if search_key == 'uniprot' or len(query.split('_')) == 2:
        u = list(UniProt.objects.filter(uniprot_identifier=query.strip()))
        if len(u) < 1:
            response = HttpResponse()
            response.write("UniProt identifier '%s' not found.\n" % query)
            return response
        else:
            return HttpResponseRedirect('/phylofacts/sequence/UniProt/%s/'% query.upper())

def citation(request):
    return render_to_response('/phylofacts/citations.html')

def dummy_search(request):
    """Dummy Search Stub

    The accession search had some commented code from ancient Limahuli.  To get
    rid of that mess and keep the accession_search view clean, yet preserve
    logic that may otherwise be lost, I've created this stub now.
    """

    keyword_search = "search_key=%s" % search_key

    r1 = Family.objects.filter(
            sequence_headers__sequence_header__uniprot__id__exact=search_key)
    r2 = Family.objects.filter(
            sequence_headers__sequence_header__header__contains=search_key)

    results = r1 | r2

    query_description = "Keyword search on '%s'" %  search_key

    return render_to_response('phylofacts/search_results.html', {
        'title': 'Search Page',
        'results': results,
    })

    # Temp move during migration/lots of moving around
    return FamilyDataGrid(request, queryset=joined_results, title="Results",
               extra_context={'extra_query': keyword_search}).\
                   render_to_response('keyword_search/results.html',
                       {'query_description': query_description})

def make_aligned_seqs_same_length(aligned_query, aligned_hit):
    realigned_query = array('c')
    realigned_hit = array('c')
    i = 0
    j = 0
    def get_to_aligned_region(i, j):
        while i < len(aligned_query) and j < len(aligned_hit) and \
              (aligned_query[i] != '-' and not aligned_query[i].isupper() \
              or aligned_hit[j] != '-' and not aligned_hit[j].isupper()):
            if aligned_query[i] == '.' and aligned_query[j] == '.':
                i += 1
                j += 1
            elif (aligned_query[i] == '.' or aligned_query[i].islower()) \
                and (aligned_hit[j] == '.' or aligned_hit[j].islower()):
                realigned_query.fromstring(aligned_query[i])
                realigned_hit.fromstring(aligned_hit[j])
                i += 1
                j += 1
            elif (aligned_query[i] == '-' or aligned_query[i].isupper()):
                if aligned_hit[j].islower():
                    realigned_query.fromstring('.')
                    realigned_hit.fromstring(aligned_hit[j])
                j += 1
            elif (aligned_hit[j] == '-' or aligned_hit[j].isupper()):
                if aligned_query[i].islower():
                    realigned_query.fromstring(aligned_query[i])
                    realigned_hit.fromstring('.')
                i += 1
        return (i, j)
    def go_through_aligned_region(i, j):
        while i < len(aligned_query) and j < len(aligned_hit) and \
            (aligned_query[i] == '-' or aligned_query[i].isupper()) and \
            (aligned_hit[j] == '-' or aligned_hit[j].isupper()):
            if (aligned_query[i] != '-' or aligned_hit[j] != '-'):
                realigned_query.fromstring(aligned_query[i])
                realigned_hit.fromstring(aligned_hit[j])
            i += 1
            j += 1
        return (i, j)

    (i, j) = get_to_aligned_region(i, j)
    while (i < len(aligned_query) and j < len(aligned_hit)):
        (i, j) = go_through_aligned_region(i, j)
        (i, j) = get_to_aligned_region(i, j)
    while (i < len(aligned_query)):
        if aligned_query[i].islower():
            realigned_query.fromstring(aligned_query[i])
            realigned_hit.fromstring('.')
        i += 1
    while (j < len(aligned_hit)):
        if aligned_hit[j].islower():
            realigned_hit.fromstring(aligned_hit[j])
            realigned_hit.fromstring('.')
    return (realigned_query.tostring(), realigned_hit.tostring())

def alignment_hmm_pairwise(request):
    error_message = None
    t = loader.get_template('phylofacts/alignment_hmm_pairwise.html')
    query_title = unquote_plus(get_parameter_from_GET('qdesc', request))
    query_length = get_parameter_from_GET('qlen', request)
    identifier = get_parameter_from_GET('identifier', request)
    basename = get_parameter_from_GET('basename', request)
    family_type_str = "ghg"
    output_file = os.path.join(PHOG_BLAST_DIR, basename + '_%s_output.pkl' % family_type_str)
    print output_file
    f = open(output_file)
    family_match_dict = cPickle.load(f)
    f.close()
    try:
        family_id = int(identifier)
    except ValueError:
        error_message = 'Invalid identifier'
        c = RequestContext(request, { 'error_message':error_message })
        return HttpResponse(t.render(c))
    family = Family.objects.get(id = family_id)
    def find_length_regions(seq):
        dedotted_seq = str(seq).translate(trivial_translation, dotdash)
        seq_length = len(dedotted_seq)
        seq_regions = contiguous_lowercase.findall(dedotted_seq)
        if dedotted_seq[0].isupper():
            seq_region = '1-'
        else:
            seq_region = '%d-' % len(seq_regions[0])
        if dedotted_seq[-1].isupper():
            seq_region = seq_region + str(len(dedotted_seq))
        else:
            seq_region = seq_region + str(len(dedotted_seq)
                                                  - len(seq_regions[-1]))
        return (seq_length, seq_region)
    print family_match_dict.keys(), family_id
    aligned_query = family_match_dict[family_id].aligned_hit
    consensus_seq = TreeNodeConsensus.objects.get(
                        tree_node =
                        family.canonical_root_node()).sequence.chars
    hmm = HMM.objects.get(tree_node = family.canonical_root_node(),
                          hmm_type='HMMER3')
    aligned_subject = consensus_seq[family_match_dict[family_id].hmm_from-1:
                                    family_match_dict[family_id].hmm_to]
    (aligned_query, aligned_subject) \
        = make_aligned_seqs_same_length(aligned_query, aligned_subject)
    subject_length = hmm.length
    identifier = family.get_accession()
    header = family.canonical_root_node().get_description(ortholog_type =
    OrthologTypes.PHOG_T_Custom, threshold=10000.0)
    subject_title = identifier
    subject_description = header
    subject_length = len(consensus_seq)
    alignment_length = len(aligned_query)
    class_of_column = {}
    identities = 0
    positives = 0
    query_specs = []
    subject_specs = []
    for i in xrange(len(aligned_query)):
        query_spec = {}
        subject_spec = {}
        query_spec['residue'] = aligned_query[i]
        subject_spec['residue'] = aligned_subject[i]
        if aligned_query[i].isupper() and aligned_subject[i].isupper():
            score = blosum62_of_residues(aligned_query[i], aligned_subject[i])
            if score >= 3:
                query_spec['class'] = 'align_high'
                subject_spec['class'] = 'align_high'
            elif score >= 1.5:
                query_spec['class'] = 'align_moderate'
                subject_spec['class'] = 'align_moderate'
            elif score >= 0.5:
                query_spec['class'] = 'align_low'
                subject_spec['class'] = 'align_low'
            if score > 0.0:
                positives += 1
            if aligned_query[i] == aligned_subject[i]:
                identities += 1
        query_specs = query_specs + [query_spec]
        subject_specs = subject_specs + [subject_spec]
    c = RequestContext(request, {
        'alignment_title': 'User Query vs %s' % subject_title,
        'query_title': query_title,
        'query_length': query_length,
        'subject_title': subject_title,
        'subject_description': subject_description,
        'subject_length': subject_length,
        'e_value': family_match_dict[family_id].i_evalue,
        'score': family_match_dict[family_id].bit_score,
        'query_region': '%s-%s' % (family_match_dict[family_id].seq_from,
                                  family_match_dict[family_id].seq_to),
        'subject_region': '%s-%s' % (family_match_dict[family_id].hmm_from,
                                  family_match_dict[family_id].hmm_to),
        'identities': identities,
        'positives': positives,
        'alignment_length': alignment_length,
        'query_specs': query_specs,
        'subject_specs': subject_specs,
        })
    return HttpResponse(t.render(c))

def sequence_search_results(request, id):
    if os.path.exists(os.path.join('/clusterfs/ohana/software/webserver/submitted_jobs/',id)):
        return render_to_response('phylofacts/sequence_search_results.html',{'id' : id})
    else:
        raise Http404

page_ze_re = re.compile('"pageSize":(\d+)')
page_num_re = re.compile('"pageNum":(\d+)')
def alignment_ajax(request, tree_node_id):
    pageSize = 20
    pageNum = 1
    if request.is_ajax():
        if request.method == 'POST':
            if u'_gt_json' in request.POST:
                m1 = page_size_re.search(request.POST[u'_gt_json'])
                m2 = page_num_re.search(request.POST[u'_gt_json'])
                if m1:
                    pageSize=int(m1.group(1))
                if m2:
                    pageNum=int(m2.group(1))
    response = HttpResponse(mimetype='application/json')
    response['Content-Disposition'] = 'attachment; filename=aln.json'
    response.write(get_alignment_rows(TreeNode.objects.get(id = int(tree_node_id)),
                                      pageSize=pageSize,pageNum=pageNum))
    return response

def serve_pdb(request, pdb_id):
    file_path = PDB.objects.get(id=pdb_id).file_path()
    if not os.path.exists(file_path):
        logging.error('File %s does not exist' % file_path)
        raise Http404
    f = open(file_path)
    response = HttpResponse(mimetype="text/plain")
    response['Content-Disposition'] = 'attachment; filename=%s.pdb.gz' \
               % (pdb_id, )
    response.write(f.read())
    f.close()
    return response

def pdb_jmol(request, pdb_id):
    return render_to_response('phylofacts/structure.html', {'pdb_id': pdb_id})

def serve_family_file_by_extension(family_accession, extension,
                                  attachment=True, qualifier = None):
    if not os.path.exists(get_dir_of_family_accession(family_accession)):
        logging.error('Directory of %s does not exist' % family_accession)
        raise Http404
    if qualifier is None:
        file_path = os.path.join(get_dir_of_family_accession(family_accession),
                               '%s.%s' % (family_accession,  extension))
    else:
        file_path = os.path.join(get_dir_of_family_accession(family_accession),
                             '%s_%s.%s' % (family_accession, qualifier, extension))
    if not os.path.exists(file_path):
        logging.error('File %s does not exist' % file_path)
        raise Http404
    f = open(file_path)
    response = HttpResponse(mimetype="text/plain")
    if attachment:
        if qualifier is None:
            response['Content-Disposition'] = 'attachment; filename=%s.%s' \
              % (family_accession, extension)
        else:
            response['Content-Disposition'] = 'attachment; filename=%s_%s.%s' \
              % (family_accession, qualifier, extension)
    response.write(f.read())
    f.close()
    return response

def family_summary_tree_json(request, family_accession, threshold_str,
                              method='ml'):
    try:
        threshold = float(threshold_str)
    except ValueError:
        raise Htttp404
    response = HttpResponse(mimetype='application/json')
    response['Content-Disposition'] = 'inline; filename=%s_%s_summary_T%s.json' \
        % (family_accession, method, threshold_str)
    response.write(annotatedSummaryTree(family_accession,
                    tree_method=method, threshold=threshold))
    return response

def ml_summary_tree_json(request, family_accession, threshold):
    return family_summary_tree_json(request, family_accession, threshold,
                                      method='ml')

def nj_summary_tree_json(request, family_accession, threshold):
    return family_summary_tree_json(request, family_accession, threshold,
                                      method='nj')

def seed(request, family_accession):
    response = HttpResponse(mimetype="text/plain")
    family = Family.objects.get(id = int(family_accession[4:]))
    response.write('>%s\n' % family.seed_sequence_header.header)
    response.write(family.seed_sequence_header.sequence.chars)
    return response


def consensus(request, family_accession):
    if not os.path.exists(get_dir_of_family_accession(family_accession)):
        logging.error('Directory of %s does not exist' % family_accession)
        raise Http404
    file_path = os.path.join(get_dir_of_family_accession(family_accession),
                             '%s.con.fa' % (family_accession))
    if not os.path.exists(file_path):
        logging.error('File %s does not exist' % file_path)
        raise Http404
    f = open(file_path)
    response = HttpResponse(mimetype="text/plain")
    # skip the header line
    f.readline()
    response.write('>lcl|%s_consensus\n' % family_accession)
    response.write(f.read())
    return response

def alignment(request, family_accession):
    return serve_family_file_by_extension(family_accession, 'a2m')

def sam_hmm(request, family_accession):
    return serve_family_file_by_extension(family_accession, 'mod')

def hmmer3_hmm(request, family_accession):
    return serve_family_file_by_extension(family_accession, 'hmm')

def family_sequences(request, family_accession):
    return serve_family_file_by_extension(family_accession, 'seqs.fa')

def family_phog_data(request, family_accession):
    return serve_family_file_by_extension(family_accession, 'phogs.txt')

def family_summary_data(request, family_accession):
    return serve_family_file_by_extension(family_accession, extension='dat',
                                          qualifier='summary')

def tree(request, family_accession, method=None):
    try:
        family = Family.objects.get(id = int(family_accession[4:]))
    except Family.DoesNotExist:
        raise Http404
    if method:
        try:
            tree = Tree.objects.get(family = family, method = method)
        except Tree.DoesNotExist:
            raise Http404
    else:
        try:
            tree = Tree.objects.get(family = family, method = 'ml')
            method = 'ml'
        except Tree.DoesNotExist:
            # There's no ML tree, try the NJ tree
            try:
                tree = Tree.objects.get(family = family, method = 'nj')
                method = 'nj'
            except Tree.DoesNotExist:
                raise Http404
    filename = "%s.%s" % (family_accession, method)
    response = HttpResponse(mimetype="text/plain")
    response['Content-Disposition'] = 'attachment; filename=%s' % filename
    try:
        tree.write_newick_to_handle(response)
    except IOException:
        raise Http404

    return response

def nj_tree(request, family_accession):
    return tree(request, family_accession, method='nj')

def ml_tree(request, family_accession):
    return tree(request, family_accession, method='ml')

def archaeopteryx_config(request):
    # This view is here so we can eventually generate this file programmatically.
    # But for now we will redirect to the static file.
    response = HttpResponse(mimetype="text/plain")
    response['Content-Disposition'] = 'attachment; filename=archaeopteryx_config.txt'
    file = open('/clusterfs/ohana/software/webserver/production/pfacts_static/static/java/archaeopteryx/archaeopteryx_config.txt','r')
    response.write(file.read())
    return response

def pfam(request, pfam_accession):
    try:
        web_dict = open('/clusterfs/ohana/bpg/pfam_pf_families/%s/web_dict'
                        % pfam_accession).read()
        web_dict = web_dict.strip()
        web_dict = eval(web_dict)
        web_dict['num_bpgs'] = len(web_dict['family_info']['values'])        
    except:
        return render_to_response('404.html',
                                  {"error_message" : "The PhyloFacts-Pfam page for %s is unavailable at this time." % pfam_accession})
    return render_to_response('phylofacts/pfam.html', { 'web_dict' : web_dict})

def biocyc_reaction(request, biocyc_reaction_id):
    try:
        web_dict = open('/clusterfs/ohana/bpg/biocyc_rxn_pf_families/%s/phog_web_dict'
                        % biocyc_reaction_id).read()
        web_dict = web_dict.strip()
        web_dict = eval(web_dict)
        web_dict['biocyc_reaction_id'] = biocyc_reaction_id
        # Getting the biocyc db for reaction
        biocyc_reactions = MetacycReaction.objects.filter(rxn_id=biocyc_reaction_id)
        preferred_dbs = ['ecoli', 'human']
        for reaction in biocyc_reactions:
            if reaction.db in preferred_dbs:
                web_dict['reaction_obj'] = reaction
                break
        web_dict['reaction_obj'] = reaction
        # Since the following information is not stored in phog_web_dict, read web_dict
        # to get the number of families.
        family_web_dict = open('/clusterfs/ohana/bpg/biocyc_rxn_pf_families/%s/web_dict'
                        % biocyc_reaction_id).read()
        family_web_dict = eval(family_web_dict)
        web_dict['num_bpgs'] = len(family_web_dict['family_info']['values'])
        # Also since we now want to add show/hide on the uniprot descriptions, modify the web_dict
        # so that the descriptions are actually a list
        for idx, val in enumerate(web_dict['phog_info']['values']):
            val['uniprot_descriptions'] = val['uniprot_descriptions'].split(';')
            web_dict['phog_info']['values'][idx] = val
    except:
        return render_to_response('404.html',
                                  {"error_message" : "The PhyloFacts-BioCyc reaction page for %s is unavailable at this time." % biocyc_reaction_id})
    return render_to_response('phylofacts/biocyc_rxn.html', { 'web_dict' : web_dict})

def sequence_data(request, family_accession):
    # This function has the beginning of a rewrie above as new_sequence_data.
    # TODO finish the rewrite above.
    family_id = int(family_accession[4:])
    # this takes care of improper accessions and 'bad' families    
    try:
        family = Family.find_by_family_id(family_id)
    except Family.DoesNotExist:
        raise Http404

    sequence_data = {}
    sequence_data['family_accession'] = family_accession
    sequence_data['final_results'] = family.get_member_data()
    return render_to_response('phylofacts/sequence_data.html', sequence_data)

def sequence_query(request):
    return render_to_response('phylofacts/sequence_query.html')

def treenode_view(request, treenode):
    ''' This populates the treenode template '''

    try:
        treenode_object = TreeNode.objects.get(id = treenode)
    except TreeNode.DoesNotExist:
        raise Http404

    # Does this node belong to an active, not bad family.

    if ((treenode_object.tree.family.status == 'bad') or not (treenode_object.tree.family.active)):
        raise Http404

    if treenode_object.is_leaf():
        if treenode_object.sequence_header.uniprot:
            return redirect('/phylofacts/sequence/UniProt/%s/' % treenode_object.sequence_header.uniprot.accession)
        else:
            raise Http404

    family_type = treenode_object.tree.family.family_type_id
    family_accession = treenode_object.tree.family.get_accession()
    node_id = "%07d_%05d" % (treenode_object.tree.id, treenode_object.left_id)
    
    # Load the species stuff
    
    # Load the MSA
    msa = treenode_object.uc_msa()
    
    # Load the domain architectures
    domain_architectures = treenode_object.get_pfam_domain_architecture_in_members()
    colormap = color_domain_architectures(domain_architectures)        

    # Load the descriptions    
    description_array = []

    # ortholog dictionary for the orthologs tab
    ortholog_info = treenode_object.orthology_supported_ancestors()
    # Load BioCyc
    biocyc_info = treenode_object.get_biocyc_reactions_and_pathways()
    go_info = treenode_object.get_go_data_with_uniprot_members()
    num_GO_annotations = go_info['total_go_annotations']
    return render_to_response('phylofacts/treenode.html', { 'treenode':treenode_object,
                                                            'page_title': 'Title', 
                                                           'msa':msa,
                                                           'taxon':treenode_object.get_taxon(),    
                                                           'num_species': treenode_object.get_num_contained_species(),
                                                           'num_structures': 0,
                                                           'num_architectures': len(domain_architectures),
                                                           'architectures' : domain_architectures,
                                                           'biocyc_info' : biocyc_info,
                                                           'domain_colormap': colormap,
                                                           'num_papers': 0,
                                                            'ortholog_info' : ortholog_info,
                                                            'go_info' : go_info,
                                                           'num_GO_annotations': num_GO_annotations,
                                                           'num_user_annotations':0},
                                                        context_instance=RequestContext(request))

def makana_admin(request):
    if request.user.is_authenticated() and request.user.is_staff:
        return render_to_response('phylofacts/fatcat_jobs_log.html', {}, context_instance=RequestContext(request))
    else:
        return render_to_response('phylofacts/login.html', {}, context_instance=RequestContext(request))

def ohana_all_jobs_table(request):
    if request.user.is_authenticated() and request.user.is_staff:
        oj = OhanaJobs.objects.all().order_by('-job_create_time')
        oj = list(oj)[:1000]
        return render_to_response('phylofacts/ohana_all_jobs_table.html', {'ohana_jobs':oj}, context_instance=RequestContext(request)) 
    else:
        raise Http404

def phylofacts_password_reset(request):
    return render_to_response('phylofacts/password_reset.html', {}, context_instance=RequestContext(request))

def my_phylofacts_view(request):
    if request.user.is_authenticated():
        fatcat_jobs =  list(FatcatJob.objects.filter(submitter_user_account=request.user))
        return render_to_response('phylofacts/my_phylofacts.html', {'fatcat_jobs': fatcat_jobs}, context_instance=RequestContext(request))
    else:
        return render_to_response('phylofacts/login.html', {}, context_instance=RequestContext(request))

def phylofacts_login(request):
    if 'member_id' in request.session:
        pass
    return render_to_response('phylofacts/login.html', {})

def phylofacts_create_user(request):
    return render_to_response('phylofacts/user_create.html',{})

def queue_status(request):
    # loads the queue status page.
    if request.user.is_staff:
        return render_to_response('phylofacts/queue_status.html', {}, context_instance=RequestContext(request))
    else:
        return render_to_response('phylofacts/login.html', {}, context_instance=RequestContext(request))

def queue_table(request):
    # table for queue status
    s = OhanaQueueStatus.objects.get(id=1)
    ret = eval(s.status)
    total_jobs = 0
    total_running = 0
    total_queued = 0
    for queue_dict in ret.values():
        for job_dict in queue_dict.values():
            if job_dict['job_status'] == 'Running':
                total_running += 1
            elif job_dict['job_status'] == 'Queued':
                total_queued += 1
            total_jobs += 1
    return render_to_response('phylofacts/queue-table.html', {'queue': ret,'total_jobs':total_jobs, 'running_jobs':total_running, 'queued_jobs':total_queued})

def fatcat_orthologs(request, id):
    try:
        job = FatcatJob.objects.get(id = id)
    except:
        raise Http404
    filepath = '/clusterfs/ohana/software/webserver/temp/fatcat/%d/' % int(id)
    filename = 'Job_%d_orthologs.csv' % int(id)
    file = open(filepath + filename, 'r')
    response = HttpResponse(mimetype="text/plain")
    response['Content-Disposition'] = 'attachment; filename=%s' % (filename)
    response.write(file.read())
    file.close()
    return response
    

def fatcat_family_orthoscope(request, family_id):
    try:
        job_family = FatcatJobFamily.objects.get(id = family_id)
    except:
        raise Http404

    if job_family.informative:
        subtree_node = job_family.informative
    else:
        subtree_node = job_family.best_tree_node
    return render_to_response('phyloscope/phyloscope_fatcat.html', 
        {
           'json_tree': job_family.family.canonical_root_node().phyloscope_json(),
           'family_accession': job_family.best_tree_node.tree.family.get_accession(),
           'hide_subfam_controls': False,
           'highlight_left_ids': "13 40", # % (job_family.best_tree_node.left_id),
           #'subtree_left_id': "%s" % (subtree_node.left_id),
           'color_left_ids': "",
           'tree_method': 'ml'
        })

def new_phyloscope(request,id):
    try:
        t = TreeNode.objects.get(id=id)
    except:
        raise Http404

    return render_to_response('phylofacts/new_phylogram2.html', {
                               'id': id,
                               'root_accession': t.get_treenode_names()['treenode_accession'],
                               'root_description': t.get_treenode_names()['treenode_name'],
                               'root_string': "%d sequences from %s" % (len(t.get_included_leaves()), t.get_taxon().scientific_name)
})

def treenode_orthoscope(request, treenode):
    try:
        treenode_object = TreeNode.objects.get(id = treenode)
    except TreeNode.DoesNotExist:
        raise Http404

    return render_to_response('phyloscope/phyloscope.html', 
            {
                'json_tree': treenode_object.phyloscope_json(),
                'family_accession': treenode_object.tree.family.get_accession(),
                'hide_subfam_controls': True,
                'tree_method': 'ml',
                'treenode': treenode_object.id
            })


def fatcat_treenode_orthoscope(request, id):
    try:
        fatcat_family = FatcatJobFamily.objects.get(id=id)
        enclosing_clade = fatcat_family.best_tree_node.get_informative_node(fatcat_family=fatcat_family)
        top_scoring_node = fatcat_family.best_tree_node
    except:
        raise Http404
    if top_scoring_node.is_leaf():
        return render_to_response('phyloscope/phyloscope2.html', 
            {
                'json_tree': enclosing_clade.phyloscope_json(top_scoring_node=top_scoring_node.id),
                'family_accession': enclosing_clade.tree.family.get_accession(),
                'hide_subfam_controls': True,
                'tree_method': 'ml',
                'top_scoring_node_id': top_scoring_node.sequence_header.id
            })
    else:
        return render_to_response('phyloscope/phyloscope2.html', 
            {
                'json_tree': enclosing_clade.phyloscope_json(top_scoring_node=top_scoring_node.id),
                'family_accession': enclosing_clade.tree.family.get_accession(),
                'hide_subfam_controls': True,
                'tree_method': 'ml',
                'top_scoring_node_id': 'TSN'
            })
   
def treenode_species_tree_json(request, treenode):
    try:
        treenode_object = TreeNode.objects.get(id = treenode)
    except TreeNode.DoesNotExist:
        raise Http404

    # Does this node belong to an active, not bad family.

    if ((treenode_object.tree.family.status == 'bad') or not (treenode_object.tree.family.active)):
        raise Http404
    
    return HttpResponse(treenode_object.get_species_tree_json(), mimetype="application/json")

def treenode_full_alignment(request, treenode):
    try:
        treenode_object = TreeNode.objects.get(id = treenode)
    except TreeNode.DoesNotExist:
        return Http404

    # Does this node belong to an active, not bad family.

    if ((treenode_object.tree.family.status == 'bad') or not (treenode_object.tree.family.active)):
        return Http404

    response = HttpResponse(mimetype="text/plain")

    response['Content-Disposition'] = 'attachment; filename=NODE%d_%d.a2m' % (treenode_object.tree.id, treenode_object.left_id)

    response.write(treenode_object.msa())

    return response

def treenode_summary_alignment(request, treenode, threshold, seed_sequence_header_id):
    """ This will send the summary alignment as an attachment.  This can
        be merged with the full alignment view at some time. """
    
    try:
        treenode_object = TreeNode.objects.get(id = treenode)
    except TreeNode.DoesNotExist:
        return Http404

    # Does this node belong to an active, not bad family.

    if ((treenode_object.tree.family.status == 'bad') or not (treenode_object.tree.family.active)):
        return Http404

    seed_treenode_queryset = TreeNode.objects.filter(tree=treenode_object.tree, sequence_header__id = seed_sequence_header_id)

    seed_treenode_object = seed_treenode_queryset[0]

    ancestor_list = seed_treenode_object.ancestors_starting_from_node(treenode_object)

    chosen_treenode = None

    for ancestor in ancestor_list:
        if (ancestor.minimum_pairwise_identity_belvu >= (float(threshold))):
            chosen_treenode = ancestor
            break

    if not chosen_treenode:
        chosen_treenode = treenode_object

    response = HttpResponse(mimetype="text/plain")

    response['Content-Disposition'] = 'attachment; filename=summary_msa_t%s_NODE%d_%d' % (threshold, treenode_object.tree.id, treenode_object.left_id)

    response.write(chosen_treenode.msa())

    return response
    

def color_domain_architectures(domain_architectures):
    colors = ('#08d9ec', '#8cf504', '#fe0000', '#06f213', '#0a4de8', '#097be9', '#09a9ea', '#08d9ec', 
             '#08edd1', '#07eea2', '#07ef73', '#06f143', '#29f305', '#5af405', '#8cf504', 
             '#bef703', '#f9cd02', '#fa9b02', '#fc6801', '#fd3401') 

    color_index = 0
    colormap = {}
    # add protein length into the dictionary, and make a dictionary of every domain and its associated color.  This allows us get all equilavent pfam domains to be the same color across all of the images.
    # this could be more clever.
    for domain_architecture in domain_architectures:
        for member in domain_architecture[1]:
            member['protein_length'] = UniProt.find_by_accession_or_identifier(member['uniprot_accession']).sequence_length    
            for domain in member['pfam_coordinates']:
                if domain[0] not in colormap:
                    colormap[domain[0]] = colors[color_index % (len(colors) - 1)]
                    color_index += 1
    return colormap

def family(request, family_accession = None):
    family_object = Family.find_by_family_id(int(family_accession.lower().replace('bpg', "")))
    treenode = family_object.canonical_root_node().id
    return redirect('/phylofacts/tree_node_view/' + str(treenode))
    """if "family_accession" not in request.GET and family_accession is None:
        return render_to_response('phylofacts/family.html',
            {'form' : FamilyAccessionForm()})
    else:
        if family_accession is None:
            family_accession = request.GET["family_accession"].strip()
            return redirect('/phylofacts/family/' + family_accession)
        form = FamilyAccessionForm({"family_accession" : family_accession})
        if not form.is_valid():
            # Error message(s) auto-generated by Django; no need for 'messages'
            # key in dict sent with response.
            return render_to_response('phylofacts/family.html',
                {'form' : form})
        else:
            try:
                id = int(family_accession[3:])
            except ValueError:
                return render_to_response('phylofacts/family.html',
                    {'form' : FamilyAccessionForm(),
                     'messages' : [messages.default_error_message]})
            try:
                family = Family.find_by_family_id(id)
            except Family.DoesNotExist:
                if "family_accession" in request.GET:
                    return render_to_response('phylofacts/family.html',
                        {'form' : form,
                         'messages' : [messages.family_does_not_exist]})
                else:
                    raise Http404

    # If we get to this point, it means 1) there were no validation errors and
    # 2) the family exists and 3) the family is not 'bad'.

    # testing family page phog table loading

    summary_data = {'family' : family, 'form' : form }
    summary_data.update({'default_num_go_items_to_display' : 10})
    summary_data.update(family.get_summary_data())
#    summary_data.update({})



    return render_to_response('phylofacts/family.html', summary_data)
    """
def family_phog_table(request, family_accession):

    try:
        id = int(family_accession[3:])
    except ValueError:
        return render_to_response('phylofacts/family.html',
            {'form' : FamilyAccessionForm(),
             'messages' : [messages.default_error_message]})
    try:
        family = Family.find_by_family_id(id)
    except Family.DoesNotExist:
        if "family_accession" in request.GET:
            return render_to_response('phylofacts/family.html',
                {'form' : form,
                 'messages' : [messages.family_does_not_exist]})
        else:
            raise Http404

    return render_to_response('phog/table.html', {
        'phogs': family.get_superorthologous_phog_rows()
    })

def uniprot_orthology_groups(request, acc_or_ident):
    try:
        object = UniProt.find_by_accession_or_identifier(acc_or_ident)
    except:
        return render_to_response('phylofacts/new_uniprot_detail_orthologs.html',
                                { 'messages': 'Your sequence is not included in any Orthology Groups.' })
    orthology_groups = object.get_orthology_groups()
    if orthology_groups:
        return render_to_response('phylofacts/new_uniprot_detail_orthologs.html',
                                { 'orthology_groups': orthology_groups,
                                  'uniprot_identifier': object.uniprot_identifier })
    else:
        return render_to_response('phylofacts/new_uniprot_detail_orthologs.html',
                                { 'messages': 'Your sequence is not included in any Orthology Groups.' })

def uniprot(request, acc_or_ident = None):
    if "acc_or_ident" not in request.GET and acc_or_ident is None:
        return render_to_response('phylofacts/uniprot_detail.html',
            {'form' : UniProtForm()}, context_instance=RequestContext(request))
    else:
        if acc_or_ident is None:
            acc_or_ident = request.GET["acc_or_ident"].strip().upper()
        form = UniProtForm({'acc_or_ident' : acc_or_ident})
        if not form.is_valid():
            # Error message(s) auto-generated by Django; no need for 'messages'
            # key in dict sent with response.
            return render_to_response('phylofacts/uniprot_detail.html',
            {'form' : form})
        else:
            # If we get to this point, the query is in a valid format. This
            # does not mean that  a UniProt object will be found in our db.
            try:
                object = UniProt.find_by_accession_or_identifier(acc_or_ident)
                ghg_families = object.ghg_families()
                pfam_family_map = object.pfam_families()
                metacyc_obj = MetacycNode.get_membership_dict(acc_or_ident)
                biocyc_obj = {'reactions': [],
                            'pathways': []}    
                if object.pfam_architecture():
                    dom_colormap = color_domain_architectures([(object.pfam_architecture()[0], [object.pfam_architecture()[1]])])
                else:
                    dom_colormap = {}
            except UniProt.DoesNotExist:
                if "acc_or_ident" in request.GET:
                    return render_to_response('phylofacts/uniprot_detail.html',
                        {'messages' : [messages.no_uniprot_record_found],
                         'form' : form,})
                else:
                    raise Http404
            return render_to_response('phylofacts/new_uniprot_detail.html',
                {'form' : form,
                 'object': object,
                 'domain_colormap': dom_colormap,
                 'pfam_family_map' : pfam_family_map,
                 'ghg_families' : ghg_families,
                 'metacyc_obj' : metacyc_obj,
                 'biocyc_info' : biocyc_obj,
                 'structures': object.get_pdb_structures(),
                 'num_user_annotations': 0,
                 'num_GO_annotations': len([a.annotation for a in summarize(object.annotations('go_biological_process'))]) +
                                    len([a.annotation for a in summarize(object.annotations('go_molecular_function'))]) +
                                    len([a.annotation for a in summarize(object.annotations('go_cellular_component'))]),
                 'go_biological_process': [a.annotation for a in summarize(object.annotations('go_biological_process'))],
                 'go_molecular_function': [a.annotation for a in summarize(object.annotations('go_molecular_function'))],
                 'go_cellular_component': [a.annotation for a in summarize(object.annotations('go_cellular_component'))],
                 }, context_instance=RequestContext(request))

def uniprot_phogs(request, acc_or_ident):
    return render_to_response('phylofacts/ajax/uniprot_phogs.html',
        {'object': get_up_by_acc_or_ident_or_404(acc_or_ident)})

def uniprot_orthologs(request, acc_or_ident, threshold=0.296874):
    obj = get_up_by_acc_or_ident_or_404(acc_or_ident)
    context = {}
    ortholog_type = PHOG.get_threshold_constants(threshold)[2]
    context_for_datagrid = {}

    ###
    newGET = request.GET.copy()
    newGET.update({'sequence_id': obj.uniprot_identifier, 'ortholog_type': ortholog_type, 'threshold': threshold})
    request.GET = newGET
    context_for_datagrid['extra_query'] = '&'.join('%s=%s' % (str(x), str(y)) for x,y in newGET.items())
    context_for_datagrid = newGET
    (ownUniProt, phogs, best_phogs,
        orthologs, phog_of_ortholog, error) = \
    getOrthologQuerySet(obj.accession, ortholog_type, threshold)
    oform = OrthologForm(request.GET)
    context['form'] = oform
    if (error != None):
        context['db_error'] = error
    elif not oform.is_valid():
        context['db_error'] = 'Error: invalid request'
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
                      for phog in best_phogs]
        context['phogs'] = phog_rows
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

        return OrthologDataGrid(request, orthologs,
              uniprot = ownUniProt,
              ortholog_type = ortholog_type,
              threshold = threshold,
              phog_of_ortholog = phog_of_ortholog,
              extra_context=context_for_datagrid).render_to_response(
                                          'phylofacts/ajax/uniprot_orthologs.html', context)

def index(request):
    operations = ['plus', 'minus', 'times']
    operation = random.sample(operations,1)[0]
    numbers = {1:'one', 2:'two', 3:'three', 4:'four',5:'five',6:'six',7:'seven',8:'eight',9:'nine',10:'ten'}
    num1 = random.randint(1,10)
    num2 = random.randint(1,10)
    question = "What is " + numbers[num1] + " " + operation + " " + numbers[num2] + "?"
    if operation == 'plus':
        answer = num1 + num2
    if operation == 'minus':
        answer = num1 - num2
    if operation == 'times':
        answer = num1 * num2
    return render_to_response('phylofacts/index_new.html', { 'captcha_question': question, 'answer_hash': hashlib.sha1(str(answer)).hexdigest()
    }, context_instance=RequestContext(request))


def viewalignment(request, familyid):
    # get alignment
    alnstr = FamilyAlignment.objects.get(family=familyid).alignment

    # replace sequence_header ids with fasta headers
    alnarr = alnstr.split('\n')
    alnstr = ''

    for line in alnarr:
        if line and line[0] == '>':
            headerid = int(line[1:len(line)])
            header   = SequenceHeader.objects.get(pk=headerid).header
            alnstr += '>' + header + '\n'

        else:
            alnstr += line + '\n'

    return render_to_response('phylofacts/alignment.html', {'alignment': alnstr})



def viewtree(request, familyid):
    # get tree
    treearr = FamilyTree.objects.filter(family=familyid)

    treestr = ''

    for treeobj in treearr:
        method = treeobj.tree_method
        tree   = treeobj.tree

        # replace sequence_header ids with fasta headers
        seqids = re.findall(r'seqh\d+', tree)

        for seqid in seqids:
            myseqid = int(seqid[4:])
            mynewid = SequenceHeader.objects.get(pk=myseqid).header
            tree    = tree.replace('%s:' % seqid, '"%s":' % mynewid)

        if method == 'ml':
            mlmodel = FamilyMlmodel.objects.get(family=familyid).model
            treestr += '[%s, model: %s] %s' % (method, mlmodel, tree)

        else:
            treestr += '[%s] %s\n' % (method, tree)

    return render_to_response('phylofacts/tree.html', {'tree': treestr})






def sequence(request, sequence_type, identifier):
    """Temporary sequence page as we are developing

    This is an incredibly simple redirect from the

    /phylofacts/sequence/UniProt/<identifier>/

    path to the:

    /phog/orthologs/<dentifier>/

    path.

    A small bit of logic was added to determine if the identifier we
    were given was an accession or a uniprot identifier.

    This code will be replaced with a proper sequence page when it is
    finished.
    """

    if sequence_type != "UniProt":
        return HttpResponse("We only support UniProt sequence pages at this time.")

    if '_' in identifier:
        return HttpResponseRedirect('/phog/orthologs/%s/' % identifier)

    else:
        u = list(UniProt.objects.filter(accession=identifier))
        if len(u) == 0:
            raise Http404
        if len(u) == 1:
            return HttpResponseRedirect('/phog/orthologs/%s/' % u[0].uniprot_identifier)
        else:
            return HttpResponse("Woah! We have %d entries for %s" % (len(u), identifier))


def download(request, accession_id, component):
    return render_to_response('phylofacts/download.html', {
        'accession_id': accession_id,
        'component': component,
    })

def downloads(request):
    return render_to_response('phylofacts/downloads.html')

def genome(request):
    return render_to_response('phylofacts/e_coli_genome.html')

def structure(request):
    return render_to_response('phylofacts/structure.html', {})

def coverage(request):
    context = {}
    context['coverage'] = {}
    context['coverage']['Eukaryotes'] = []
    context['coverage']['Bacteria'] = []
    context['coverage']['Archaea'] = []
    eukaryotes = UniProtTaxonomy.objects.get(scientific_name = 'Eukaryota')
    bacteria = UniProtTaxonomy.objects.get(scientific_name = 'Bacteria')
    archaea = UniProtTaxonomy.objects.get(scientific_name = 'Archaea')
    def is_in_taxon(ancestor, descendant):
        return (descendant.left_id >= ancestor.left_id and
                descendant.right_id <= ancestor.right_id)
    files = glob.glob('/clusterfs/ohana/bpg/coverage/redundant/pfam/after_QFO8/ID/QFO/*.coverage')
    for file in files:
        name = os.path.splitext(os.path.split(file)[1])[0]
        name_components = name.split('_')
        taxon_id = int(name_components[0])
        taxon = UniProtTaxonomy.objects.get(id = taxon_id)
        taxon_name = taxon.scientific_name
        f = open(file)
        lines = f.readlines()
        f.close()
        num_covered = int(lines[2].split()[1].strip('(').strip(')'))
        num_covered += int(lines[3].split()[1].strip('(').strip(')'))
        num_uncovered = int(lines[4].split()[1].strip('(').strip(')'))
        total_num = num_covered + num_uncovered
        percent_covered = float(num_covered) / total_num * 100.0
        coverage_str = '%0.1f%%' % percent_covered
        if num_covered > 0:
            if is_in_taxon(eukaryotes, taxon):
                coverage_list = context['coverage']['Eukaryotes']
            elif is_in_taxon(archaea, taxon):
                coverage_list = context['coverage']['Archaea']
            else:
                coverage_list = context['coverage']['Bacteria']
            coverage_list.append((taxon.left_id, taxon_name, total_num,
                                        num_covered, coverage_str))
    for domain in context['coverage'].keys():
        context['coverage'][domain].sort()
    return render_to_response('phylofacts/coverage.html', context)

def functional_site_prediction(request):
    return render_to_response('phylofacts/functional_site_prediction_index.html')

def functional_site_prediction_job(request, jobID):
    job = FxnSitePredictionJob.objects.get(id=jobID)
    fasta_sequence = wrap(job.fasta_sequence)
    return render_to_response('phylofacts/functional_site_prediction_job.html', {'job':job, 'fasta_sequence': fasta_sequence })

def functional_site_prediction_job_results(request, jobID):
    # parse the intrepid file and 
    #functional_site_context = 
    job_directory = os.path.join("/clusterfs/ohana/software/webserver/temp/", jobID)
    print job_directory
    if not os.path.exists(job_directory):
        raise Http404

    (intrepid_datatable_rows, intrepid_cons_js_score, intrepid_res_name,
        intrepid_cons_js_rank) = parse_intrepid_output(job_directory)

    return render_to_response('phylofacts/intrepid_job_results.html', 
                    {   
                        'intrepid_datatable_rows': intrepid_datatable_rows,
                        'intrepid_cons_js_score': intrepid_cons_js_score,
                        'intrepid_res_name': intrepid_res_name,
                        'intrepid_cons_js_rank': intrepid_cons_js_rank,                                             'id': jobID
                    })
# To serve the files for functional site prediciton

def serve_FSP_results(jobID, shortname):
    filename = os.path.join("/clusterfs/ohana/software/webserver/temp/", jobID, shortname)
    f = open(filename)
    response = HttpResponse(mimetype="text/plain")
    response['Content-Disposition'] = 'attachment; filename=%s' % shortname
    response.write(f.read())
    f.close()
    return response

def functional_site_prediction_summary_msa(request, jobID):
    return serve_FSP_results(jobID, "job_%d_summary.msa" % int(jobID))

def functional_site_prediction_summary_tree(request, jobID):
    return serve_FSP_results(jobID, "job_%d_summary.tree" % int(jobID))

def functional_site_prediction_full_msa(request, jobID):
    return serve_FSP_results(jobID, "job_%d.msa" % int(jobID))

def functional_site_prediction_full_tree(request, jobID):
    return serve_FSP_results(jobID, "job_%d.tree" % int(jobID))

def functional_site_prediction_hmm(request, jobID):
    return serve_FSP_results(jobID, "job_%d.hmm" % int(jobID))

def functional_site_prediction_intrepid_rank(request, jobID):
    return serve_FSP_results(jobID, "output.rank")

def functional_site_prediction_intrepid_score(request, jobID):
    return serve_FSP_results(jobID, "output.aux")

#### From old 'book' views
#"""All views related to /books/..... URLs"""

# For Operation Obese Ocelot
def fatcat(request):
    return HttpResponseRedirect('/fatcat/')

def fatcat_example_one(request):
    response_file = '/clusterfs/ohana/software/webserver/examples/fatcat/F7G2Z0/job2121.httpResponse'
    file = open(response_file, 'r')
    resp = HttpResponse(file.read())
    file.close()
    return resp

def fatcat_example_two(request):
    response_file = '/clusterfs/ohana/software/webserver/examples/fatcat/Q13530/job2108.httpResponse'
    file = open(response_file, 'r')
    resp = HttpResponse(file.read())
    file.close()
    return resp

def fatcat_job_new(request, id):
    try:
        job = FatcatJob.objects.get(id=id)
    except:
        raise Http404

    fasta_sequence = wrap(job.fasta_sequence)
    if job.status.id < 9:
        return render_to_response('phylofacts/fatcat_job_new.html', {'job':job, 'sequence':fasta_sequence}, context_instance=RequestContext(request))
    else:
        job_dir = "/clusterfs/ohana/software/webserver/temp/fatcat/" + str(job.id) + "/"
        response_file = job_dir + "job%d.httpResponse" % job.id
        if os.path.exists(response_file):
            file = open(response_file, "r")
            resp = HttpResponse(file.read())
            file.close()
            return resp
        else:
            return render_to_response('phylofacts/fatcat_job_results.html', { 'job': job, 'fasta_sequence': fasta_sequence }, context_instance=RequestContext(request))
        
def fatcat_job(request, id):
    job = FatcatJob.objects.get(id=id)
    #To display it all pretty like
    fasta_sequence = wrap(job.fasta_sequence)
    if job.status.id < 9:
        return render_to_response('phylofacts/fatcat_job.html', { 'job': job, 'fasta_sequence': fasta_sequence }, context_instance=RequestContext(request))
    else:
        #return render_to_response('phylofacts/fatcat_job_results.html', { 'job': job, 'fasta_sequence': fasta_sequence })
        job_dir = "/clusterfs/ohana/software/webserver/temp/fatcat/" + str(job.id) + "/"
        response_file = job_dir + "job%d.httpResponse" % job.id
        if os.path.exists(response_file):
            file = open(response_file, "r")
            resp = HttpResponse(file.read())
            file.close()
            return resp
        else:
            return render_to_response('phylofacts/fatcat_job_results.html', { 'job': job, 'fasta_sequence': fasta_sequence }, context_instance=RequestContext(request))

def fatcat_about_detail(request, stage_number=None):
    return render_to_response('phylofacts/fatcat_about.html', {})

def fatcat_about(request):
    return render_to_response('phylofacts/fatcat_about.html', {})

def fatcat_precalc(request):
    return render_to_response('phylofacts/fatcat_precalc.html', {})

def fatcat_supplementary(request):
    return render_to_response('phylofacts/fatcat_supplementary.html', {})

    
def fatcat_help(request):
    return render_to_response('phylofacts/fatcat_help.html', {}, context_instance=RequestContext(request))

def fatcat_faq(request):
    return render_to_response('phylofacts/fatcat_faq.html', {})

def fatcat_family(request, family_id):
    family = FatcatJobFamily.objects.get(id=family_id)
    return render_to_response('phylofacts/fatcat_family.html', { 'family': family })

def fatcat_family_tree(request, family_id):
    '''Added separate fatcat tree since the query is needed for getting an informative node,
    for functions like get pwid, coverage etc.'''
    response = HttpResponse(mimetype='application/xml')
    fatcat_family = FatcatJobFamily.objects.get(id=family_id)
    tree_node = fatcat_family.best_tree_node
    informative_node = tree_node.get_informative_node(fatcat_family=fatcat_family)
    if informative_node:
        response.write(informative_node.tree.family.canonical_root_node().phyloxml(
                annotations=['GO', 'EC', 'uniprot_accession', 'uniprot_identifier', 
                    'uniprot_description', 'PHOG', 'taxonomy', 'SFLD'],
                    stop_on=["KERF"], unless=[informative_node.id], threshold=70, 
                highlight=[tree_node.id]))
    else:
        response.write(TreeNode.objects.get(id = int(tree_node.id)).tree.family.canonical_root_node().phyloxml(
                annotations=['EC', 'uniprot_accession', 'uniprot_identifier', 
                    'uniprot_description', 'PHOG', 'SFLD', 'taxonomy'], 
                stop_on=["KERF"], threshold=70,
                highlight=[tree_node.id]))
    return response

def tree_node(request, id):
    return render_to_response('phylofacts/tree_viewer.html', { 'id': id })

def tree_node_tree(request, id):
    response = HttpResponse(mimetype='application/xml')
    tree_node = TreeNode.objects.get(id = int(id))
    response.write(TreeNode.objects.get(id = int(id)).phyloxml(
        annotations=['EC', 'consensus_descriptions', 'GO', 'swissProt', 'literature', 
                     'structure', 'percent_id', 'mrca', 'taxonomy']))
    return response
    
def tree_node_phyloxml_tree_download(request, id):
    try:
        tree_node = TreeNode.objects.get(id = int(id))
    except:
        raise Http404

    response = HttpResponse(tree_node.phyloxml(), mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.phyloxml' % tree_node.get_treenode_names()['treenode_accession']

    return response
    
def tree_node_masked_msa_download(request, id):
    try:
        tree_node = TreeNode.objects.get(id = int(id))
    except:
        raise Http404
    
    response = HttpResponse(tree_node.uc_msa(), mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.afa' % tree_node.get_treenode_names()['treenode_accession']

    return response

def tree_node_a2m_msa_download(request, id):
    try:
        tree_node = TreeNode.objects.get(id = int(id))
    except:
        raise Http404
    
    response = HttpResponse(tree_node.msa(), mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.a2m' % tree_node.get_treenode_names()['treenode_accession']

    return response

def tree_node_hmm_download(request, id):
    try:
        tree_node = TreeNode.objects.get(id = int(id))
    except:
        raise Http404
     
    response = HttpResponse(tree_node.subtree_hmm(), mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.hmm' % tree_node.get_treenode_names()['treenode_accession']
    
    return response

def tree_node_members_download(request, id): 
    try:
        tree_node = TreeNode.objects.get(id = int(id))
    except:
        raise Http404

    response = HttpResponse(tree_node.get_member_fasta(), mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.members' % tree_node.get_treenode_names()['treenode_accession']

    return response

def tree_node_newick_tree_download(request, id):
    try:
        tree_node = TreeNode.objects.get(id = int(id))
    except:
        raise Http404

    response = HttpResponse(tree_node.newick(), mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.newick' % tree_node.get_treenode_names()['treenode_accession']

    return response

# pointer to the static file
def benchmarking_pdf(request):
    #return render_to_response('phylofacts/fatcat_about.html', {})
    #response = HttpResponse(tree_node.newick(), mimetype='text/plain')
    #response['Content-Disposition'] = 'attachment; filename=%s.newick' % tree_node.get_treenode_names()['treenode_accession']

    return HttpResponseRedirect('/static/doc/benchmarking.pdf')
    
def tutorial_pdf(request):
    return HttpResponseRedirect('/static/doc/FATCAT_Tutorial.pdf')
