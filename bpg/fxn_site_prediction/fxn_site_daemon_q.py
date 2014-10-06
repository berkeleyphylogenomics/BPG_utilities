#!/clusterfs/ohana/software/bin/python2.7
"""
This is the functional site prediction job daemon.  We have to do this
to allow webusers access to our queue.
"""

import sys, re, os, time, subprocess
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
    submit_script_file = open("qsub_script.sh", "w")
    submit_script_file.write("#!/bin/sh\n\n" + 
                             ("#PBS -N Functional_Site_Prediction_Job_%s\n" % job.id) + 
                             ("#PBS -q web\n") + 
                             ("#PBS -l nodes=1:ppn=1,walltime=1000\n") +
                             ("#PBS -o %s\n" % JOB_LOG_FILE) + 
                             ("#PBS -e %s\n" % JOB_LOG_FILE) + 
                             ("#PBS -v job_id=%s\n" % job.id) + 
                             ("#PBS -r y\n\n") + 
                             ("python /clusterfs/ohana/home/cyrus_afrasiabi/ohana_repository/bpg/fxn_site_prediction/fxn_site_prediction.py\n"))
    submit_script_file.close()
    args = ["qsub", os.path.abspath(submit_script_file.name)]
    print args
    process = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    pbs_job_id = (process.communicate()[0]).split(".")[0]
    job.status_id = 2
    job.save()

    logger.info("Submitting %s to the grid to get functional site predictions with id %s" % (job.id, pbs_job_id))

    if pbs_job_id:
        job.pbs_job_id = pbs_job_id
        job.save()
    
    os.remove(os.path.abspath(submit_script_file.name))
 
    return pbs_job_id

while (1 == 1):
    logger.info("-- Polling")
    for job in FxnSitePredictionJob.objects.filter(status__id=1):
        submit_fxn_site_prediction_job(job)

    # We check every 10 seconds
    time.sleep(10)
