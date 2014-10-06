#!/clusterfs/ohana/software/bin/python2.7
#PBS -z

import os, sys,time
sys.path.append('/home/jdobbie/ohana_repository')
sys.path.append('/clusterfs/ohana/software/lib/python2.7/site-packages/')
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
from pfacts003.phylofacts.models import FatcatJob
from pfacts003.utils.hmm import top_k_families, top_family_per_region

job_id = os.environ["job_id"]
job = FatcatJob.objects.get(id=job_id, status__id=2)

if not job:
    raise

job.status_id = 3
job.save()

top_family_hmm_rows = top_k_families(job.fasta, 1000000, type="ghg")
for hmm_row in top_family_hmm_rows:
    id = hmm_row.target_name[4:]
    fatcat_family = job.families.create(
        family_id=id,
        family_e_value = hmm_row.full_e_value,
        family_type = 'Multiple Domain Architecture')

    fatcat_family.coords.create(
            ali_from=hmm_row.ali_from,
            ali_to=hmm_row.ali_to
        )

#    for dom_row in hmm_row.doms:
#        fatcat_family.coords.create(
#            ali_from=dom_row.ali_from,
#            ali_to=dom_row.ali_to
#        )
#        break # only want the first one I guess
        

top_family_hmm_rows = top_k_families(job.fasta, 1000000, type="pfam")
for hmm_row in top_family_hmm_rows:
    id = hmm_row.target_name[4:]
    fatcat_family = job.families.create(
        family_id=id,
        family_e_value = hmm_row.full_e_value,
        family_type = 'Pfam')

    fatcat_family.coords.create(
            ali_from=hmm_row.ali_from,
            ali_to=hmm_row.ali_to
        )

#    for dom_row in hmm_row.doms:
#        fatcat_family.coords.create(
#            ali_from=dom_row.ali_from,
#            ali_to=dom_row.ali_to
#        )
#        break # only want the first one I guess
#seen = set()
#for hmm_row in top_family_per_region(job.fasta, type="pfam"):
#for hmm_row in top_k_families(job.fasta, 1000000, type="pfam"):
#    id = int(hmm_row.target_name[4:])
#
#    if id not in seen:
#        seen.add(id)
#        family = job.families.create(
#            family_id=id,
#            family_e_value = hmm_row.full_e_value,
#            family_type = 'Pfam')
#    else:
#        family = job.families.get(family=id)
#
#    for dom_row in hmm_row.doms:
#        fatcat_family.coords.create(
#            ali_from=dom_row.ali_from,
#            ali_to=dom_row.ali_to
#        )
#        break # only want the first one I guess
##    family.coords.create(
##        ali_from=hmm_row.ali_from,
##        ali_to=hmm_row.ali_to
##    )
    


job.status_id = 4
job.save()
