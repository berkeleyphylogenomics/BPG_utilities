#!/clusterfs/ohana/software/fatcat/bin/python
import os
import sys

os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
sys.path.append('/clusterfs/ohana/software/fatcat/lib/python2.7/site-packages/')
#sys.path.insert(0,'/home/cyrus_afrasiabi/ohana_repository')

import lockfile, logging, time
from daemon import DaemonContext

import pbs, smtplib
from pfacts003.fatcat.models import FatcatJob
from email.MIMEText import MIMEText
from email.mime.multipart import MIMEMultipart

JOB_LOG_FILE = '/clusterfs/ohana/software/fatcat/logs/fatcat2_job_log'
THIS_EMAIL = 'phylofacts.webmaster@gmail.com'

def submit_fatcat_job(job):
    # This is how we are passing the fasta and job id to the script
    server_name = pbs.pbs_default()
    c = pbs.pbs_connect(server_name)
    attropl = pbs.new_attropl(8)

    attropl[0].name  = pbs.ATTR_N
    attropl[0].value = "FAT-CAT Job: %s" % job.id

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

    attropl[6].name = pbs.ATTR_m
    attropl[6].value = 'abe'

    attropl[7].name = pbs.ATTR_M
    attropl[7].value = THIS_EMAIL

    job_id = pbs.pbs_submit(c, attropl, "/clusterfs/ohana/software/fatcat/scripts/fatcat2_pipeline.py", 'web', 'NULL')
    logger.info("Submitting %s to the grid with id %s" % (job.id, job_id))

    if job_id: 
        job.pbs_job_id = job_id
        job.status_id = 3
        job.save()
    else:
        pass

    pbs.pbs_disconnect(c)

    return job_id

with DaemonContext():
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler('/clusterfs/ohana/software/fatcat/logs/fatcat2_daemon_log'))
    logger.setLevel(logging.INFO)
    logger.info("-- Starting Daemon Loop")
    while (1 == 1):
        logger.info("-- Polling")
        for job in FatcatJob.objects.filter(status__id=1, test=False):
            submit_fatcat_job(job)
        time.sleep(10)
