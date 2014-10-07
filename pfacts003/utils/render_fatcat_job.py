#!/usr/bin/env python



import os, sys
sys.path.insert(0,'/home/cyrus_afrasiabi/ohana_repository')
from pfacts003.phylofacts.models import FatcatJob
from django.shortcuts import render_to_response
from textwrap import wrap

os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"

def main():
    if len(sys.argv) != 6:
        print "Wrong number of arguments, correct usage ./render_fatcat_job.py [fatcat job id] [mda coverage (float 0-1)] [pfam coverage (float 0-1)] [min query to node %ID (integer 1-99)] [min pairwise %ID of subtree (integer 1-99)]"
        return 1
    fatcat_job_id = int(sys.argv[1])
    mda_coverage = float(sys.argv[2])
    pfam_coverage = float(sys.argv[3])
    min_query_to_node_pid = float(sys.argv[4])
    min_pairwise_pid = float(sys.argv[5])
    
    try:
        fatcat_job = FatcatJob.objects.get(id = fatcat_job_id)
        fatcat_job.criteria_get_best_nodes_mda_overlap = mda_coverage
        fatcat_job.criteria_get_best_nodes_pfam_overlap = pfam_coverage
        fatcat_job.criteria_get_best_nodes_subtree_min_pid = min_pairwise_pid
        fatcat_job.criteria_get_best_nodes_query_to_subtree_hmm_pid = min_query_to_node_pid
        fatcat_job.save()
    except:
        print "Fatcat job %d doesn't exist" % fatcat_job_id
        return 1
    job_dir = "/clusterfs/ohana/software/webserver/temp/fatcat/%d/" % fatcat_job_id
    job_html_file = job_dir + "job%d.httpResponse" % fatcat_job_id
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)
    if os.path.exists(job_html_file):
        os.remove(job_html_file)
    file = open(job_html_file, "w")
    file.write(render_to_response('phylofacts/fatcat_job_results.html',
                {'job': fatcat_job, 'fasta_sequence': wrap(fatcat_job.fasta_sequence)}).content)
    file.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
