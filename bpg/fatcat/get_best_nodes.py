#!/clusterfs/ohana/software/bin/python2.7
#PBS -z

import os,sys,time
sys.path.append('/home/jdobbie/ohana_repository')
sys.path.append('/clusterfs/ohana/software/lib/python2.7/site-packages/')
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
from pfacts003.phylofacts.models import FatcatJob, FatcatJobFamily
from pfacts003.utils.hmm import hmmscan, parse_tblout

job_id = os.environ["job_id"]
job = FatcatJob.objects.get(id=job_id, status__id=5)

if not job:
    raise

job.status_id = 6
job.save()

for fatcat_family in FatcatJobFamily.objects.filter(fatcat_job=job):
    hmm_file = fatcat_family.family.subtree_hmm_file
    (normal, dom) = hmmscan(job.fasta, hmm_file, [ '--domE', '1e-4', '-E', '1e-3' ])
    normal_table = parse_tblout(normal)
    normal_table.sort(reverse=True)
    if len(normal_table) > 0:
        fatcat_family.best_tree_node_id = normal_table[0].target_name
        fatcat_family.best_node_e_value = normal_table[0].full_e_value
        fatcat_family.save()

job.status_id = 7
job.save()
