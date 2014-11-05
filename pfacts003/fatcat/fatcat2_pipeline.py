#!/clusterfs/ohana/software/bin/python2.7
#PBS -z

import os, sys, time
sys.path.append('/clusterfs/ohana/software/prod/lib/python2.4/site-packages/')
sys.path.append('/clusterfs/ohana/software/lib/python2.7/site-packages/')
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
#sys.path.insert(0,'/home/cyrus_afrasiabi/ohana_repository')
#sys.path.insert(0,'/home/ddineen/ohana_repository')
from pfacts003.fatcat.models import FatcatJob, FatcatJobFamily, FatcatJobFamilyMatchData
from ete2 import Tree

job_id = os.environ['job_id']
job = FatcatJob.objects.get(id=job_id)


job.status_id = 3
job.save()

#starting scan vs pfam hmms
job.scan_pfams()

job.scan_ghmms(type='Multiple Domain Architecture')

job.status_id = 4
job.save()

job.scan_ghmms(type='Pfam')

if not job.get_families_passing_stage1():
    job.status_id = 12
    job.log_job_completion()
    sys.exit()

stage_2_families = []

if job.run_ellis_island:
    pfacts_families = job.pass_through_ellis_island()
#    pfacts_families = job.taxon_based_ellis_island()

    if not pfacts_families:
        job.status_id = 12
        job.log_job_completion()
        sys.exit()
    stage_2_families = [FatcatJobFamily.objects.get(fatcat_job = job, family = fam) for fam in pfacts_families]
else:
    stage_2_families = job.families.all()

job.status_id = 5
job.save()

job.scan_subtree_hmms_mt(stage_2_families)

job.status_id = 6
job.save()

job.analyze_top_scoring_nodes()

job.status_id = 7
job.save()

msa = job.build_fatcat_job_msa()

if msa is None:
    job.status_id = 12
    job.log_job_completion()
    sys.exit()

job.status_id = 8
job.save()

# build the tree
#if job.has_enclosing:
full_tree = job.build_fatcat_job_tree(msa)
    
try:
    # root the tree here using ete
    tree_object = Tree(full_tree)
    outgroup = tree_object.get_midpoint_outgroup()
    tree_object.set_outgroup(outgroup)

    # cut the tree
    [sub_tree, sub_msa] = job.cut_fatcat_job_tree(msa, tree_object.write())
    # populate the members
    job.status_id = 9
    job.save()
        
    job.populate_fatcat_job_members(sub_msa, sub_tree)
except:
    job.status_id = 9
    job.save()
    job.populate_fatcat_job_members(msa, full_tree)

    
job.status_id = 10
job.save()

job.analyze_job_members_for_orthology()

job.status_id = 11
job.save()

#job.store_functions_transferred(gene=job.predict_query_gene_name(), description=job.predict_query_description())
job.store_functions_transferred(gene=job.predict_query_gene_name(), description=job.predict_query_description(), ec=job.predict_query_ec_number())

job.populate_member_go_annotations()

job.status_id = 12
job.save()

job.log_job_completion()
