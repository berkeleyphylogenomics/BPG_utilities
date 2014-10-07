import sys,os, gzip
from django.shortcuts import render_to_response
from django.http import Http404, HttpResponse
from django.template import RequestContext
from pfacts003.intrepid.models import IntrepidJob, IntrepidStatus, IntrepidJobPfamDomain
from pfacts003.intrepid.consts import *
from pfacts003.intrepid.utils import repairPDB

def download(request, id, type):
    """ handles the file download """
    def get_domain(did, jid):
        # returns the domain specified by did, should have the same job id as jid
        try:
            d = IntrepidJobPfamDomain.objects.get(id=did)
        except:
            return None

        if d.intrepid_job.id != jid:
            return None
        return d

    def get_response(p, fn):
        # build the download file
        f = open(p, 'rb')
        response = HttpResponse(f.read(), mimetype='text/plain')
        response['Content-Disposition'] = 'attachment; filename=%s' % fn
        f.close()
        return response
    
    # parse the type field
    try:
        job = IntrepidJob.objects.get(id=int(id))
        (download_type, domid) = type.split('_')
        domid = int(domid)
    except:
        raise Http404
        
    if (download_type == 'tree'):
        if (domid) == 0:
            response = get_response(job.tree_path, 'FSP_Job%d_tree.newick' % job.id)
        else:
            d = get_domain(domid, job.id)
            if d:
                response = get_response(d.tree_path, 'FSP_Job%d_%s_tree.newick' % (job.id, d.pfam_shortname))
            else:
                raise Http404
    elif (download_type == 'msa'):
        if (domid) == 0:
            # full length msa
            response = get_response(job.homolog_alignment_tree_path, 'FSP_Job%d_member_alignment.afa' % job.id)
        else:
            d = get_domain(domid, job.id)
            if d:
                # domain msa
                response = get_response(d.homolog_alignment_tree_path, 'FSP_Job%d_%s_member_alignment.afa' % (job.id, d.pfam_shortname))
            else:
                raise Http404
    elif (download_type == 'scores'):
        if (domid) == 0:
            # full length scores
            response = get_response(job.intrepid_aux_file_path, 'FSP_Job%d_INTREPID_Residue_Scores.txt' % job.id)
        else:
            d = get_domain(domid, job.id)
            if d:
                # domain scores
                response = get_response(d.intrepid_aux_file_path, 'FSP_Job%d_%s_INTREPID_Residue_Scores.txt' % (job.id, d.pfam_shortname))
            else:
                raise Http404
    elif (download_type == 'discern'):
        response = get_response(job.discern_output_file_path, 'FSP_Job%d_DISCERN_Scores.txt' % (job.id))
    elif (download_type == 'progout'):
        response = get_response(job.program_outputs_path, 'FSP_Job%d_program_outputs.txt' % (job.id))
    else:
        raise Http404
    return response

def serve_this_job_pdb(request, id):
    try:
        job = IntrepidJob.objects.get(id=id)
    except:
        raise Http404
    
    response = HttpResponse(mimetype='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s.pdb' % (job.name)
    f = open(job.original_pdb_file_path, 'rb')
    if not job.pdb_id[4:]:
        response.write(repairPDB(f.read(), options=['-chain', 'A']))
    else:
        response.write(f.read())
    f.close()
    return response    

def serve_pdb(request, pid):
    # serves the appropriate pdb file, if there is an underscore in the pid, it tries to serve the zip file.
    if "_" in pid:
        type = pid.split('_')[0]
        pid = pid.split('_')[1]
    pdbid = pid[:4]
    chainid = pid[4:]
    if not os.path.exists(os.path.join(PDB_LIBRARY_PATH, pdbid[1:3], 'pdb%s.ent.gz' % (pdbid.lower()))):
        raise Http404
    response = HttpResponse(mimetype='text/plain')
    if type == "Z":
        f = open(os.path.join(PDB_LIBRARY_PATH, pdbid[1:3], 'pdb%s.ent.gz' % (pdbid.lower())), 'rb')
        response['Content-Disposition'] = 'attachment; filename=%s.pdb.gz' % (pdbid.lower())
    else:
        f = gzip.open(os.path.join(PDB_LIBRARY_PATH, pdbid[1:3], 'pdb%s.ent.gz' % (pdbid.lower())), 'rb')
        response['Content-Disposition'] = 'attachment; filename=%s.pdb' % (pdbid.lower())
    response.write(f.read())
    f.close()
    return response

def index(request):
    """handles displaying the initial protein analysis form,
    the progress of running jobs
    """ 
    return render_to_response('intrepid/index.html', {
            'uniprot_all_val': UNIPROT_ALL,
            'uniref100_val': UNIREF100,
            'uniref90_val': UNIREF90,
            'blast_new_val': BLAST_NEW_VERSION,
            'blast_old_val': BLAST_OLD_VERSION
        }, context_instance = RequestContext(request))

def results(request, id):
    """pulls up a display of existing results"""
    try:
        job = IntrepidJob.objects.get(id=int(id))
    except:
        raise Http404
    
    # Should we redirect to the progress page?
    if not job.is_done:
        return render_to_response('intrepid/progress.html', {
                    'job': job,
                    'progress_update_interval': PROGRESS_PAGE_TIMEOUT_INTERVAL
                }, context_instance = RequestContext(request))
    else:
        if job.is_discern_job:
            return render_to_response('intrepid/discern_results.html', {
                    'job':job,
                    'query_length': len(job.fasta_sequence),
                    'substatus': job.get_substatus_text()
                }, context_instance = RequestContext(request))
        else:
            return render_to_response('intrepid/intrepid_results.html', {
                    'job':job,
                    'query_length': len(job.fasta_sequence),
                    'substatus': job.get_substatus_text()
                }, context_instance = RequestContext(request))
