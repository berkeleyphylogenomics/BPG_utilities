#!/clusterfs/ohana/software/bin/python2.7
"""
This is the functional site prediction job daemon.  We have to do this
to allow webusers access to our queue.
"""

import sys, re, os, time
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
sys.path.append('/home/cyrus_afrasiabi/ohana_repository')
sys.path.append('/clusterfs/ohana/software/lib/python2.7/site-packages/')

import pbs 
from pfacts003.phylofacts.models import FxnSitePredictionJob

import logging
logger = logging.getLogger()
logger.addHandler(logging.FileHandler('/dev/stdout'))
logger.setLevel(logging.INFO)

JOB_LOG_FILE = '/home/cyrus_afrasiabi/ohana_repository/bpg/fxn_site_prediction/log'

def submit_fxn_site_prediction_job(job):
    # This is how we are passing the fasta and job id to the script
    server_name = pbs.pbs_default()
    c = pbs.pbs_connect(server_name)

    print server_name
    print c
    
    attropl = pbs.new_attropl(7)

    attropl[0].name  = pbs.ATTR_N
    attropl[0].value = "Functional Site Prediction Job: %s" % job.id

    attropl[1].name  = pbs.ATTR_l
    attropl[1].resource = 'nodes'
    attropl[1].value = '1:ppn=1'

    attropl[2].name  = pbs.ATTR_o
    attropl[2].value = JOB_LOG_FILE

    attropl[3].name  = pbs.ATTR_e
    attropl[3].value = JOB_LOG_FILE

    attropl[4].name  = pbs.ATTR_v
    attropl[4].value = "job_id=%s" % (job.id)

    attropl[5].name  = pbs.ATTR_r
    attropl[5].value = 'y'

    attropl[6].name  = pbs.ATTR_l
    attropl[6].resource = 'walltime'
    attropl[6].value = '1000'

    job.status_id = 2
    job.save()

    job_id = pbs.pbs_submit(c, attropl, "/home/cyrus_afrasiabi/ohana_repository/bpg/fxn_site_prediction.py", 'web', 'NULL')
    logger.info("Submitting %s to the grid to get functional site predictions with id %s" % (job.id, job_id))

    if job_id:
        job.pbs_job_id = job_id
        job.save()

    pbs.pbs_disconnect(c)
    
    return job_id

while (1 == 1):
    logger.info("-- Polling")
    for job in FxnSitePredictionJob.objects.filter(status__id=1):
        submit_fxn_site_prediction_job(job)

    # We check every 10 seconds
    time.sleep(10)
