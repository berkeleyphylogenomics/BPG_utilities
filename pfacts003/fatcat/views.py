# Create your views here.
from django.http import Http404, HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.db.models import Q
from pfacts003.fatcat.models import FatcatJob, FatcatJobFamily, FatcatJobMembers
from textwrap import wrap
from pfacts003.fatcat.consts import *
from StringIO import StringIO
from Bio import Phylo

def index(request):
    return render_to_response('fatcat/index.html', {}, context_instance = RequestContext(request)) 

def results(request, id):
    try:
        job = FatcatJob.objects.get(id=id)
    except:
        raise Http404

    fasta_sequence = wrap(job.fasta_sequence)
    if not job.is_done:
        return render_to_response('fatcat/progress.html', {'job': job, 'num_fatcat_stages':NUM_FATCAT_STAGES, 'sequence': fasta_sequence}, context_instance = RequestContext(request))
    else:
        # this can be done better....
        taxon = set()
        for member in job.get_orthologs():
            taxon.add(member.uniprot.taxon.id)
        try:
            functions = eval(job.functions)
        except:
            descriptions = ''
            genes = ''
            ec = ''
        else:
            try:
                descriptions = functions['description']
            except:
                descriptions = ''
            try:
                genes = functions['gene']
            except:
                genes = ''
            try:
                ec = functions['ec']
            except:
                ec = ''
        return render_to_response('fatcat/results.html', {
            'job':job,
            'num_stage1_families': job.families.all().count(),
            'num_stage2_families': len(job.get_families_passing_stage2()),
            'num_genomes_with_orthologs': len(taxon),
            'num_distant_clades': job.families.filter(passed_stage1=True, passed_stage2=False).count(),
            'tree_size': job.members.all().count(),
            'sequence': fasta_sequence,
            'descriptions': descriptions,
            'genes': genes,
            'ec': ec,
            'high_threshold': CONSENSUS_UNIPROT_THRESHOLD_FOR_HIGH_CONFIDENCE,
            'medium_threshold': CONSENSUS_UNIPROT_THRESHOLD_FOR_MEDIUM_CONFIDENCE
        }, context_instance = RequestContext(request))

def family_tree(request, id):
    try:
        fatcat_family = FatcatJobFamily.objects.get(id=id)
        enclosing_clade = fatcat_family.enclosing_clade_root
        top_scoring_node = fatcat_family.top_scoring_node
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

def orthologs_download(request, id):
    try:
        job = FatcatJob.objects.get(id=id)
    except:
        raise Http404

    response = HttpResponse(mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=FAT_CAT_Job_%d_orthologs.csv' % job.id
    response.write(job.get_ortholog_file(delimiter=','))
    return response
    
def tree(request, id, type='newick'):
    # get the tree for a fatcat job
    try:
        job = FatcatJob.objects.get(id=id)
    except:
        raise Http404

    if type == 'newick':
        tree = Phylo.read(StringIO(job.tree_string),'newick')
        for leaf in tree.find_clades():
            if leaf.name:
                if str('query') in leaf.name:
                    if(job.fasta_header):
                        leaf.name = '"QUERY_' + job.fasta_header[1:].replace(':', '.').replace(';', '') + '"'
                    else:
                        leaf.name = '"QUERY"'
                else:
                    for m in job.members.filter(is_query = False):
                        if m.uniprot.accession in leaf.name:
                            leaf.name = '"' + m.uniprot.uniprot_identifier + ' ' + m.uniprot.description.replace(':', '').replace(';', '') + '"'
                leaf.name = leaf.name.replace(':', '_').replace(';', '_').replace(',', '').replace(' ', '_').replace('/', '_')
        job.tree_string = tree.format("newick")

        response = HttpResponse(mimetype='text/plain')
        response['Content-Disposition'] = 'attachment; filename=FAT_CAT_Job_%d_tree.newick' % job.id
        response.write(job.tree_string)
    elif type == 'phyloxml':
        response = HttpResponse(mimetype='application/xml')
        response['Content-Disposition'] = 'attachment; filename=FAT_CAT_Job_%d_tree.phyloxml' % job.id
        response.write(job.phyloxml_string)
    else:
        raise Http404
    return response 

def alignment(request, id):
    try:
        job = FatcatJob.objects.get(id=id)
    except:
        raise Http404

    response = HttpResponse(mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=FAT_CAT_Job_%d_alignment.fasta' % job.id
    response.write(job.get_fatcat_job_msa())
    return response 
