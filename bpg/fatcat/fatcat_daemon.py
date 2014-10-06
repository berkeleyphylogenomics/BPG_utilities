#!/clusterfs/ohana/software/fatcat/bin/python
''' Greetings sparklepony!

You may be wondering why this file exists, and you may even hate me for creating
it.  Or you've read enough of my other code to hate me in general, and just 
wish for some explanation for my madness.  Well, shod off.

The purpose of this file is to poll the database for jobs that are able to run,
and to submit them to the queue.  SCS doesn't want the apache user to have this
power.  They say it is for security, but their reasons are kind of wishy washy.
Okay, they are stupid, but, we don't want to piss off SCS any more than is
required, so in this situation, the path of least resistant was to create a 
daemon that will run as me, and handle all of the queue magic.  Not that not
submitting things with apache is that bad of a thing, but they also are fighting
us on a daemon user, and lots of other stuff.
'''
import os
import sys

os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
sys.path.append('/clusterfs/ohana/software/fatcat/lib/python2.7/site-packages/')

import lockfile, logging, time
from daemon import DaemonContext

import pbs, smtplib
from pfacts003.phylofacts.models import FatcatJob
from email.MIMEText import MIMEText
from email.mime.multipart import MIMEMultipart

JOB_LOG_FILE = '/clusterfs/ohana/software/fatcat/logs/job_log'

def submit_get_families_job(job):
    # This is how we are passing the fasta and job id to the script
    server_name = pbs.pbs_default()
    c = pbs.pbs_connect(server_name)
    attropl = pbs.new_attropl(7)

    attropl[0].name  = pbs.ATTR_N
    attropl[0].value = "FAT-CAT Get Families: %s" % job.id

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

    job_id = pbs.pbs_submit(c, attropl, "/clusterfs/ohana/software/fatcat/scripts/get_families.py", 'web', 'NULL')
    logger.info("Submitting %s to the grid to get best families with id %s" % (job.id, job_id))

    if job_id: 
        job.get_families_pbs_job_id = job_id
        job.save()
    else:
        job.status_id = 1
        job.save()

    pbs.pbs_disconnect(c)

    return job_id

def submit_get_subfamilies_job(job):
    # This is how we are passing the fasta and job id to the script
    server_name = pbs.pbs_default()
    c = pbs.pbs_connect(server_name)
    attropl = pbs.new_attropl(5)

    attropl[0].name  = pbs.ATTR_N
    attropl[0].value = "FAT-CAT Get Sub-Families: %s" % job.id

    attropl[1].name  = pbs.ATTR_l
    attropl[1].resource = 'nodes'
    attropl[1].value = '1:ppn=1'

    attropl[2].name  = pbs.ATTR_o
    attropl[2].value = JOB_LOG_FILE

    attropl[3].name  = pbs.ATTR_e
    attropl[3].value = JOB_LOG_FILE

    attropl[4].name  = pbs.ATTR_v
    attropl[4].value = "job_id=%s" % (job.id)

    job.status_id = 5
    job.save()

    job_id = pbs.pbs_submit(c, attropl, "/clusterfs/ohana/software/fatcat/scripts/get_best_nodes.py", 'web', 'NULL')
    logger.info("Submitting %s to the grid to get best nodes with id %s" % (job.id, job_id))

    if job_id: 
        job.get_best_nodes_pbs_job_id = job_id
        job.save()

    pbs.pbs_disconnect(c)

    return job_id


def send_email(job):
    if job.user_email:
        EMAIL_FROM = "phylofacts.webmaster@gmail.com"
        EMAIL_PASSWORD = "lanikai324c"

        text_body = '''PhyloFacts FAT-CAT Job %d has completed.
            \n\nhttp://makana-test.berkeley.edu/phylofacts/fatcat/%s
            \n%s\n%s
            \n\n''' % (job.id, job.id, job.fasta_header, job.fasta_sequence[:10])
        html_body = '''<h2>PhyloFacts FAT-CAT Job %d has completed</h2>.
            \n\n<p><a href="http://makana-test.berkeley.edu/phylofacts/fatcat/%s">click here to view the results</a></p>
            \n<p>%s</p>\n<p>%s</p>
            \n\n''' % (job.id, job.id, job.fasta_header, job.fasta_sequence[:10])
        text = MIMEText(text_body, 'plain')
        html = MIMEText(html_body, 'html')

        msg = MIMEMultipart('alternative')
        msg['Subject'] = "PhyloFacts FAT-CAT Job %d is complete." % job.id
        msg['From'] = EMAIL_FROM
        msg['To'] = job.user_email
        msg.attach(text)
        msg.attach(html)

        # Sending email
        server = smtplib.SMTP('smtp.gmail.com:587')
        server.ehlo()
        server.starttls()
        server.ehlo()
        server.login(EMAIL_FROM, EMAIL_PASSWORD)
        server.sendmail(EMAIL_FROM, job.user_email, msg.as_string())
        server.close()

    job.status_id = 8
    job.save()

with DaemonContext():
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler('/clusterfs/ohana/software/fatcat/logs/daemon_log'))
    logger.setLevel(logging.INFO)
    logger.info("-- Starting Daemon Loop")
    while (1 == 1):
        logger.info("-- Polling")
        for job in FatcatJob.objects.filter(status__id=1):
            submit_get_families_job(job)
    
        for job in FatcatJob.objects.filter(status__id=4):
            submit_get_subfamilies_job(job)
    
        for job in FatcatJob.objects.filter(status__id=7):
            send_email(job)
    
        time.sleep(10)
