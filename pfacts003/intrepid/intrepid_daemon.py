#!/clusterfs/ohana/software/bin/python2.7
import os
import sys

os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
sys.path.append('/clusterfs/ohana/software/lib/python2.7/site-packages/')
#sys.path.insert(0,'/home/cyrus_afrasiabi/ohana_repository')

import lockfile, logging, time
from daemon import DaemonContext

import pbs, smtplib
from pfacts003.intrepid.models import IntrepidJob
from pfacts003.intrepid.consts import *
from email.MIMEText import MIMEText
from email.mime.multipart import MIMEMultipart

JOB_LOG_FILE = '/clusterfs/ohana/software/intrepid/logs/intrepid_job_log'

def submit_intrepid_job(job):
    # This is how we are passing the fasta and job id to the script
    server_name = pbs.pbs_default()
    c = pbs.pbs_connect(server_name)
    attropl = pbs.new_attropl(6)

    attropl[0].name  = pbs.ATTR_N
    attropl[0].value = "INTREPID Job: %s" % job.id

    attropl[1].name  = pbs.ATTR_l
    attropl[1].resource = 'nodes'
    attropl[1].value = '1:ppn=8'

    attropl[2].name  = pbs.ATTR_o
    attropl[2].value = JOB_LOG_FILE

    attropl[3].name  = pbs.ATTR_e
    attropl[3].value = JOB_LOG_FILE

    attropl[4].name  = pbs.ATTR_v
    attropl[4].value = "job_id=%s" % (job.id)

    attropl[5].name = pbs.ATTR_l
    attropl[5].resource = 'walltime'
    attropl[5].value = '48:00:00'
    
    if job.development_job:
        job_id = pbs.pbs_submit(c, attropl, "/clusterfs/ohana/software/intrepid/scripts/intrepid_development_pipeline.py", 'web', 'NULL')
    else:
        job_id = pbs.pbs_submit(c, attropl, "/clusterfs/ohana/software/intrepid/scripts/intrepid_pipeline.py", 'web', 'NULL')
    logger.info("Submitting %s to the grid with id %s" % (job.id, job_id))

    if job_id: 
        job.pbs_job_id = job_id
        job.status_id = JOB_SUBMITTED
        job.save()
    else:
        pass

    pbs.pbs_disconnect(c)

    return job_id

with DaemonContext():
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler('/clusterfs/ohana/software/intrepid/logs/intrepid_daemon_log'))
    logger.setLevel(logging.INFO)
    logger.info("-- Starting Daemon Loop")
    while (1 == 1):
        logger.info("-- Polling")
        for job in IntrepidJob.objects.filter(status__id=WAITING_FOR_SUBMISSION):
            submit_intrepid_job(job)
        time.sleep(10)
