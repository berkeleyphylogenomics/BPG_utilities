#!/clusterfs/ohana/software/bin/python2.7
#PBS -z

''' This pipeline is for development.  It will get submitted to the queue if the 
development_job boolean field in the intrepid_job table is True.  This is so that
we can develop and have a separate pipeline that is run using the production repository.
'''

import os, sys, time
sys.path.append('/clusterfs/ohana/software/prod/lib/python2.4/site-packages/')
sys.path.append('/clusterfs/ohana/software/lib/python2.7/site-packages/')
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
sys.path.insert(0,'/home/cyrus_afrasiabi/ohana_repository')
from pfacts003.intrepid.models import IntrepidJob 
from pfacts003.intrepid.consts import *

job_id = os.environ['job_id']
os.environ['DATA_DIR'] = INTREPID_DATA_DIR

job = IntrepidJob.objects.get(id=job_id)

job.status_id = SCAN_PFAMS
job.save()

job.scan_pfams()

job.status_id = GATHER_HOMOLOGS
job.save()

job.gather_homologs()

job.status_id = PROCESS_MSA 
job.save()

job.mask_msa()

job.status_id = BUILD_TREE
job.save()

job.build_tree()

job.status_id = RUN_INTREPID
job.save()

job.run_intrepid()

job.status_id = RUN_DISCERN
job.save()

if job.is_discern_job:
    job.run_discern()
    job.populate_discern_scores()

job.log_job_completion()

sys.exit()
