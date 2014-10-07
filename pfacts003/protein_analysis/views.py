from django.shortcuts import render_to_response
from django.http import Http404, HttpResponse
from django.template import RequestContext
from pfacts003.protein_analysis.models import ProteinAnalysisJob
from pfacts003.protein_analysis.consts import *

def index(request):
    """handles displaying the initial protein analysis form,
    the progress of running jobs
    """ 
    return render_to_response('protein_analysis/index.html', {}, context_instance = RequestContext(request))

def scanning_domains(request, id):
    ''' Shows a progress bar when we are scanning for domains '''
    try:
        job = ProteinAnalysisJob.objects.get(id=int(id))
    except:
        raise Http404

    return render_to_response('protein_analysis/scanning_domains.html', {}, context_instance = RequestContext(request))

def results(request, id):
    """pulls up a display of existing results"""
    try:
        job = ProteinAnalysisJob.objects.get(id=int(id))
    except:
        raise Http404
    
    # Should we redirect to the scanning domains page?
    if job.status < SCAN_FOR_DOMAINS_STAGES:
        return render_to_response('protein_analysis/scanning_domains.html', {
                    'job': job,
                    'progress_update_interval': PROGRESS_PAGE_TIMEOUT_INTERVAL
                }, context_instance = RequestContext(request))

    return render_to_response('protein_analysis/results.html', {
                'job':job,
                'mda_domain_id': job.domains.get(type='mda').id,
                'query_length': len(job.fasta_sequence)
            }, context_instance = RequestContext(request))
