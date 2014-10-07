from django.db import models
from django.db import connection
from django.db.models import Q
from django.contrib.auth.models import User
from pfacts003.phylofacts.models import TreeNode, Family, UniProt, UniProtTaxonomy, UniProtGene, OhanaJobs, PhyloFactsEmail, Annotation, EC, OntologySubsetMembership, UniProtEC
from pfacts003.utils.hmm import pfamAscan, hmmscan, parse_domtblout, hmmsearch, hmmbuild
from pfacts003.utils.fatcat import *
from pfacts003.utils.koa_trim_tree import *
from pfacts003.utils.alignment import mafft
from bpg.common.utils.sequenceutils.fasta_utils import get_uniprot_accession_from_header
from StringIO import StringIO
from pfacts003.fatcat.consts import *
from datetime import datetime
from Bio import AlignIO, Phylo
import tempfile, time, os, glob
import NewQueue as Queue
import Queue
import threading
#import shutil
import sys
from itertools import groupby
import pprint

try:
    import json
except:
    import simplejson as json

class FatcatJobStatus(models.Model):
    id = models.IntegerField(primary_key=True)
    status = models.TextField()

    class Meta:
        db_table = u'fatcat2_job_status'

class FatcatJob(models.Model):
    id = models.AutoField(primary_key=True)
    # FASTA 
    fasta_header = models.TextField()
    fasta_sequence = models.TextField()
    # alignment of query - only the row of the MSA corresponding to the query
    query_alignment = models.TextField()
    # timestamps
    created_at = models.DateTimeField(default=datetime.now)
    # job status
    status = models.ForeignKey(FatcatJobStatus, related_name="jobs")
    # email address and email subject
    user_email = models.TextField()
    email_subject = models.TextField()
    # job comments
    fatcat_job_name = models.TextField()
    # queue job id
    pbs_job_id = models.TextField()
    # ohana job id
    ohana_job = models.ForeignKey(OhanaJobs, related_name='fatcat2_job')
    # stores a string of the parameterization the user chose for this job
    parameterization = models.TextField()
    # job parameters
    stage1_e_value = models.FloatField()
    stage1_mda_family_query_coverage = models.FloatField()
    stage1_mda_family_hmm_coverage = models.FloatField()
    stage1_pfam_family_hmm_coverage = models.FloatField()
    stage2_e_value = models.FloatField()
    stage2_mda_family_query_coverage = models.FloatField()
    stage2_mda_family_hmm_coverage = models.FloatField()
    stage2_pfam_family_hmm_coverage = models.FloatField()
    stage2_query_hmm_pid = models.FloatField()
    stage2_orthology_methods = models.TextField()
    stage2_enclosing_clade_kerf_threshold = models.FloatField()
    stage2_redundant_methods_required = models.IntegerField()
    stage2_max_sequences = models.IntegerField()
    stage3_orthology_methods = models.TextField()
    stage3_kerf_threshold = models.IntegerField()
    stage3_query_ortholog_pid = models.FloatField()
    stage3_query_coverage = models.FloatField()
    stage3_ortholog_coverage = models.FloatField()
    stage3_ortholog_cluster_similarity = models.FloatField()
    stage3_max_sequences = models.IntegerField()
    stage4_consensus_uniprot_weighting_lambda = models.FloatField()
    stage4_consensus_uniprot_weighting_amplitude = models.FloatField()
    stage4_consensus_uniprot_threshold_for_high_confidence = models.FloatField()
    stage4_consensus_uniprot_threshold_for_medium_confidence = models.FloatField()
    # This is for storing a json object corresponding to the functions transferred to the query
    # functions
    functions = models.TextField()
    # ortholog tree
    tree_string = models.TextField()
    phyloxml_string = models.TextField()
    # should this job be saved?
    save_job = models.BooleanField()
    # seconds for job completion    
    job_completion_time = models.FloatField()
    # this isn't really necessary...
    is_done = models.BooleanField()
    # is this a test job?
    test = models.BooleanField()
    # ellis island parameters
    run_ellis_island = models.BooleanField()
    ellis_island_maxseqs = models.IntegerField()
    ellis_island_maxfams = models.IntegerField()
    ellis_island_member_coverage = models.FloatField()
    ellis_island_query_coverage = models.FloatField()
    ellis_island_taxon_clustering_pid = models.FloatField()

    @property
    def fasta(self):
        if self.fasta_header.startswith('>'):
            return self.fasta_header + '\n' + self.fasta_sequence
        else:
            return '>' + self.fasta_header + '\n' + self.fasta_sequence
  
    # These should be implemented as a function in the template that multiplies by 100
    # oh well :'( 
    @property
    def readable_stage1_mda_family_query_coverage(self):
        return int(self.stage1_mda_family_query_coverage*100)

    @property
    def readable_stage1_mda_family_hmm_coverage(self):
        return int(self.stage1_mda_family_hmm_coverage*100)

    @property
    def readable_stage1_pfam_family_hmm_coverage(self):
        return int(self.stage1_pfam_family_hmm_coverage*100)

    @property
    def readable_ellis_island_member_coverage(self):
        return int(self.ellis_island_member_coverage*100)
        
    @property
    def readable_ellis_island_query_coverage(self):
        return int(self.ellis_island_query_coverage*100)

    @property
    def readable_ellis_island_taxon_clustering_pid(self):
        return int(self.ellis_island_taxon_clustering_pid*100)        
    @property
    def readable_stage2_mda_family_query_coverage(self):
        return int(self.stage2_mda_family_query_coverage*100)

    @property
    def readable_stage2_mda_family_hmm_coverage(self):
        return int(self.stage2_mda_family_hmm_coverage*100)

    @property
    def readable_stage2_pfam_family_hmm_coverage(self):
        return int(self.stage2_pfam_family_hmm_coverage*100)
 
    ########################################################################################################
    #
    # Functions for getting various sets of families, members, etc.
    #
    ########################################################################################################

    def get_families_passing_stage1(self):
        try:
            return self._get_families_passing_stage1
        except:
            pass

        self._get_families_passing_stage1 = [x for x in self.families.all() if x.passes_stage1_coverage_conditions]
        return self._get_families_passing_stage1
    
    def get_families_passing_stage2(self):
        try:
            return self._get_families_passing_stage2
        except:
            pass

        self._get_families_passing_stage2 = [x for x in self.families.filter(passed_stage2=True)]
        return self._get_families_passing_stage2

    def get_ortholog_clusters(self):
        ''' This function returns a list of tuples, the first member being the cluster number,
            the second member being the fatcatjobmember that is the representative of that cluster. '''
        try:
            return self._get_ortholog_clusters
        except:
            pass

        self._get_ortholog_clusters = [ ( x.cluster_num, x ) for x in self.members.filter(is_cluster_rep=True, classification = 'Orthologs') ]
        return self._get_ortholog_clusters
    
    def get_other_sequence_match_clusters(self):
        ''' This function returns a list of tuples, the first member being the cluster number,
            the second member being the fatcatjobmember that is the representative of that cluster. '''
        try:
            return self._get_other_sequence_match_clusters
        except:
            pass

        self._get_other_sequence_match_clusters = [ ( x.cluster_num, x ) for x in self.members.filter(is_cluster_rep=True, classification='Other Sequence Matches') ]
        return self._get_other_sequence_match_clusters

    def get_cluster_members(self, cn):
        ''' This function will return a list of the members inside a cluster. '''
        return list(self.members.filter(cluster_num = cn))

    ########################################################################################################
    #
    # Functions that create crap in the database and should only be run once per job - to create the entries
    # maybe we should switch the creates to get_or_create....
    #
    ########################################################################################################

    def scan_pfams(self):
        pfam_hits = pfamAscan(self.fasta, options=['-E','1e-3'])
        for pfam in pfam_hits:
            fatcat_job_pfam = self.pfams.create(
                pfam_description = pfam.description,
                pfam_accession = pfam.target_accession,
                pfam_shortname = pfam.target_name,
                ali_from = pfam.ali_from,
                ali_to = pfam.ali_to,
                hmm_from = pfam.hmm_from,
                hmm_to = pfam.hmm_to,
                c_evalue = pfam.this_domain_c_value,
                i_evalue = pfam.this_domain_i_value
            )
        return None

    def scan_ghmms(self, type):
        # hmmscan vs ghg family hmms
        if type == 'Pfam':
            hmm_file = '/clusterfs/ohana/bpg/pfacts/phylofacts3_domain_090811.hmms'
        elif type == 'Multiple Domain Architecture':
            hmm_file = '/clusterfs/ohana/bpg/pfacts/phylofacts3_GHG_090811.hmms'
 
        # add options for each parameterization
        if self.parameterization == 'High recall':
            hmmscan_options = ['-E', str(self.stage1_e_value), '--domE', str(self.stage1_e_value)]
        elif self.parameterization == 'High precision':
            hmmscan_options = ['-E', str(self.stage1_e_value), '--domE', str(self.stage1_e_value)]
        elif self.parameterization == 'Remote homologs':
            hmmscan_options = ['-E', str(self.stage1_e_value), '--domE', str(self.stage1_e_value)]
        elif self.parameterization == 'Partial sequence':
            hmmscan_options = ['--max', '-E', str(self.stage1_e_value), '--domE', str(self.stage1_e_value)]
        elif self.parameterization == 'Custom':
            hmmscan_options = ['-E', str(self.stage1_e_value), '--domE', str(self.stage1_e_value)]
        else:
            hmmscan_options = ['-E', str(self.stage1_e_value), '--domE', str(self.stage1_e_value)]

        ((normal, domain), ec) = hmmscan(self.fasta, hmm_file, hmmscan_options)
        # parse the hmmscan output
        domain_hmmscan_rows = parse_domtblout(domain)
        # sort by bit score
        domain_hmmscan_rows.sort(reverse=True)
        seen = set()
        top_families = []
        # get the best families (unique) that have alignment lengths greater than the cutoff
        for row in domain_hmmscan_rows:
            if ( row.ali_to - row.ali_from ) < STAGE1_MINIMUM_ALIGNMENT_LENGTH:
                continue
            if row.target_name not in seen:
                top_families.append(row)
                seen.add(row.target_name)
        # store the families and the match data
        for fam in top_families:
            try:
                f = Family.objects.get(id=fam.target_name[4:])
            except:
                continue
                       
            if ((fam.this_domain_i_value < self.stage1_e_value) and f.canonical_root_node().get_included_leaves()):
                # get family description
                if type=='Pfam':
                    description = f.get_pfams().pop()[0].name
                else:
                    try:
                        description = f.canonical_root_node().consensus_description.get().consensus_uniprot_description
                    except:
                        description = f.canonical_root_node().get_treenode_names()['treenode_name']
                fatcat_family = self.families.create(
                    family = f,
                    family_e_value = fam.this_domain_i_value,
                    family_type = type,
                    family_description = description,
                    passed_stage1 = False,
                    passed_stage2 = False    
                )
                
                # get family taxonomy
                try:
                    taxonomy = f.canonical_root_node().mrca.get().mrca
                except:
                    taxonomy = f.canonical_root_node().get_taxon()
                
                if taxonomy is not None:
                    fatcat_family.family_taxonomy = taxonomy
                    fatcat_family.save()
                # get family mda
    
                fatcat_family.match_data.create(
                    stage1_ali_from = fam.ali_from,
                    stage1_ali_to = fam.ali_to,
                    stage1_hmm_from = fam.hmm_from,
                    stage1_hmm_to = fam.hmm_to,
                    stage1_hmm_length = fam.target_length,
                    stage1_region_aligned = int(fam.ali_to) - int(fam.ali_from),
                    stage1_query_coverage = float(fam.ali_to - fam.ali_from)/len(self.fasta_sequence),
                    stage1_hmm_coverage = float(fam.hmm_to - fam.hmm_from)/int(fam.target_length)
                )
        return None
    
    def populate_fatcat_job_members(self, msa, tree):
        # takes in an msa and tree and creates the job members
        query_length = len(self.fasta_sequence)
        query_ali = ''
        # first find the query record
        for record in msa:
            if record.id.startswith('query'):
                query_ali = str(record.seq)
                self.members.create(
                    alignment = str(record.seq),
                    is_query = True,
                    member_coverage = 0,
                    query_coverage = 0,
                    member_to_query_pid = 0,
                    is_ortholog = False,
                    reason_eliminated = ''
                )
                break
        # hopefully this should never happen
        if not query_ali:
            return None

        for record in msa:
            if not record.id.startswith('query'):
                # these are really hard to debug...just got messed by one now...
                try:
                    u = UniProt.objects.get(accession=get_uniprot_accession_from_header(record.id))
                    member_length = u.sequence_len
                    # calculate % ID and coverage
                    aligned_chars = 0
                    match = 0
                    identity = 0
                    for (index, qchar) in enumerate(query_ali):
                        mchar = record[index]
                        if (qchar.isupper() and mchar.isupper()):
                            aligned_chars += 1
                            if mchar == qchar:
                                identity += 1
                    qcov = float(100.0*aligned_chars/query_length)
                    mcov = float(100.0*aligned_chars/member_length)
                    try:
                        pwid = float(100.0*identity/aligned_chars)
                    except:
                        pwid = 0.0
                    # do this for now. it should be removed when i - or you - think more about it
                    if qcov > 100:
                        qcov = 100.0
                    if mcov > 100:
                        mcov = 100.0
                    if pwid > 100:
                        pwid = 100.0
                    self.members.create(
                        uniprot = u,
                        alignment = str(record.seq),
                        is_query = False,
                        member_coverage = mcov,
                        query_coverage = qcov,
                        member_to_query_pid = pwid,
                        is_ortholog = False, # nothing yet
                        reason_eliminated = '', # nothing yet
                        is_cluster_rep = False,
                        classification = ''
                    ) 
                except:
                    pass
        return None

    def analyze_job_members_for_orthology(self):
        members = list(self.members.filter(is_query=False))
        for member in members:
            if member.member_to_query_pid > self.stage3_query_ortholog_pid:
                if member.query_coverage > self.stage3_query_coverage:
                    if member.member_coverage > self.stage3_ortholog_coverage:
                        member.reason_eliminated = 'Passed'
                        member.classification = 'Orthologs'
                        member.save()
                    else:
                        member.reason_eliminated = 'Does not meet Stage 3 Candidate Ortholog coverage criterion'
                        member.classification = 'Other Sequence Matches'
                        member.save()
                else:
                    member.reason_eliminated = 'Does not meet Stage 3 Query coverage criterion'
                    member.classification = 'Other Sequence Matches'
                    member.save()
            else:
                member.reason_eliminated = 'Does not meet minimum Query to Candidate Ortholog % ID'
                member.classification = 'Other Sequence Matches'
                member.save()

        # first cluster the orthologs
        num_clusters, ortholog_cluster_dict = cluster_member_list(list(self.members.filter(is_query=False, classification='Orthologs').order_by('-member_to_query_pid')))
       
        seen_tax = {}
        # Now, remove all ortholog clusters in taxa where there is a higher scoring cluster...these are
        # "paralogs" - this might not be right....it might give wrong behavior....we will see.
        for tax_id, cluster in ortholog_cluster_dict.items():
            if tax_id in seen_tax:
                if cluster[0]['representative'].member_to_query_pid < seen_tax[tax_id].member_to_query_pid:
                    for m in cluster[0]['members']:
                        m.classification = 'Other Sequence Matches'
                        m.save()
                else:
                    seen_tax[tax_id] = cluster[0]['representative']
            else:
                seen_tax[tax_id] = cluster[0]['representative']
 
        # then cluster the other sequence matches
        num_clusters, osm_cluster_dict = cluster_member_list(list(self.members.filter(is_query=False, classification='Other Sequence Matches').order_by('-member_to_query_pid')), cluster_number_start=num_clusters)
        # Now go through the members in both dicts and populate their cluster numbers and representatives
        for tax_id, cluster in ortholog_cluster_dict.items():
            new_cluster = self.clusters.create(
                cluster_representative = cluster[0]['representative'],
                taxon = UniProtTaxonomy.objects.get(id=tax_id),
                cluster_size = len(cluster[0]['members']),
            )    
            for member in cluster[0]['members']:
                member.cluster = new_cluster
                if member == cluster[0]['representative']:
                    member.is_cluster_rep = True
                member.save()

        for tax_id, cluster in osm_cluster_dict.items():
            new_cluster = self.clusters.create(
                cluster_representative = cluster[0]['representative'],
                taxon = UniProtTaxonomy.objects.get(id=tax_id),
                cluster_size = len(cluster[0]['members']),
            )    
            for member in cluster[0]['members']:
                member.cluster = new_cluster
                if member == cluster[0]['representative']:
                    member.is_cluster_rep = True
                member.save()
        '''top_scoring_taxonomies = {}
        cluster_representatives = {}
        cluster_num = 1

        # this whole clustering thing should be replaced, it's quite silly.

        for member in sorted_members:
            taxon = member.uniprot.taxon.id
            if member.uniprot.taxon.id in top_scoring_taxonomies:
                if (top_scoring_taxonomies[taxon].member_to_query_pid - member.member_to_query_pid < 5) and member.uniprot.in_swissprot_f and not cluster_representatives[taxon].uniprot.in_swissprot_f:
                    cluster_representatives[taxon] = member
            else:
                top_scoring_taxonomies[taxon] = member
                cluster_representatives[taxon] = member
        for member in sorted_members:
            taxon = member.uniprot.taxon.id
            if cluster_representatives[taxon] == member:
                member.reason_eliminated = 'None'
                member.is_ortholog = True
                member.cluster_num = cluster_num 
            elif top_scoring_taxonomies[taxon].member_to_query_pid - member.member_to_query_pid < 5:
                member.reason_eliminated = 'In a %s genome cluster with %s' % (member.uniprot.taxon.scientific_name, cluster_representatives[taxon].uniprot.uniprot_identifier)
            else:
                member.reason_eliminated = 'Outside of %s cluster in genome %s with higher similarity to the query' % (cluster_representatives[taxon].uniprot.uniprot_identifier, cluster_representatives[taxon].uniprot.taxon.scientific_name)
            member.save()
        '''
        return ortholog_cluster_dict, osm_cluster_dict
    
    ######################################################################################################
    #
    # Other functions
    #
    ######################################################################################################

    def compile_fasta_for_ellis_island(self, queue):
        # this is a function that gets the chunk from the queue
        while True:
            chunk = queue.get()
            # get the fasta for these uniprot objects
            uniprot_list = chunk[0]
            this_dir = chunk[1]
            this_file = open(tempfile.mkstemp(dir=this_dir, prefix='job%d'%self.id)[1], 'w')
            count = 0
            for uniprot in uniprot_list:
                count += 1
                this_file.write(uniprot.fasta_from_db+'\n')
            # print 'fasta chunk count ' + str(count)
            this_file.close()
            connection.close()
            queue.task_done()

    def pass_through_ellis_island(self):
        ''' This implements a speed up for fatcat. 
            Returns a list of families that should be used in stage2
        '''
        TEMP_DIRECTORY = '/clusterfs/ohana/software/webserver/incoming/'
        # make a temporary directory here for this job
        this_directory = tempfile.mkdtemp(dir=TEMP_DIRECTORY)
        # this is the sequence file for the hmmscan for this job
        fasta_file = open(os.path.join(this_directory, 'ellis_island.fasta'), 'w')
        # first get the unique members in these families, and make a list of their accessions
        family_list = [fam.family.id for fam in self.get_families_passing_stage1()]
        uniprot_set = find_uniprots_for_family_list(family_list)
        # now we need to compile the FASTA for these uniprot sequences
        queue = Queue.Queue()
        for i in range(THREADS):
            t = threading.Thread(target=self.compile_fasta_for_ellis_island, args=(queue, ))
            t.setDaemon(True)
            t.start()
        for chunk in chunkify(uniprot_set, THREADS):
            if chunk:
                queue.put((chunk, this_directory), True)
        # wait for completion
        queue.join()
        temporary_fasta_paths = [file for file in glob.glob(os.path.join(this_directory,'job%d*'%self.id))]
        # concatenate these member fasta files
        args = ['cat'] + temporary_fasta_paths
        process = subprocess.Popen(args, shell=False, stdout=fasta_file, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        process.communicate()
        fasta_file_path = fasta_file.name
        fasta_file.close()
        # Now, get the top k1 sequences to the query
        # first, get an hmm of the query
        hmm_file = open(tempfile.mkstemp(dir=this_directory)[1],'w')
        hmm_file.write(hmmbuild(self.fasta, 'query')[0])
        hmm_file.close()
        ((full_length, domain), error) = hmmsearch(hmm_file.name, fasta_file_path, 
                                    options=['--cpu', '8', '-E', '13-3', '--domE','1e-3'])
        # get a list of the accessions sorted by bit score and only those that pass coverage conditions
        # this already has an implicit evalue threshold imposed by hmmer
        top_accessions = [row.target_name.split('|')[1] for row in sorted(parse_domtblout(domain), 
            key=lambda hmmerrow: hmmerrow.full_score) 
            if float(row.ali_to - row.ali_from)/len(self.fasta_sequence) >= self.ellis_island_member_coverage and
            float(row.hmm_to - row.hmm_from)/int(row.target_length) >= self.ellis_island_query_coverage]
        #print 'top accession count: ' + str(len(top_accessions))
        # get the top k1 unique accessions
        top_accessions = top_accessions[0:self.ellis_island_maxseqs]
        family_list = find_families_for_uniprot_accession_list(top_accessions, family_list)
        families_needed = []
        while top_accessions:
            this_family = ellis_island_get_best_family(top_accessions, [fam.id for fam in family_list])
            for member in get_members_in_accession_list_in_family(top_accessions, this_family):
                top_accessions.remove(member)
            family_list.remove(this_family)
            families_needed.append(this_family)
        
        # change permissions of everything in this directory, and remove it
        # not working yet (leaves empty directory)
        subprocess.call(['chmod', '-R', '777', this_directory])
        subprocess.call(['rm', '-rf', this_directory])
        #shutil.rmtree(this_directory)
        #subprocess.call(['rmdir', this_directory])

        return families_needed[0:self.ellis_island_maxfams]

    # DD: new ellis island implementation - idea is to not have multiple entries from same speceis swamp results
    # Returns a list of families that should be used in stage2
    def taxon_based_ellis_island(self):
        # make a temporary directory here for this job
        TEMP_DIRECTORY = '/clusterfs/ohana/software/webserver/incoming/'
        this_directory = tempfile.mkdtemp(dir=TEMP_DIRECTORY)
        # this is the sequence file for the hmmscan for this job
        fasta_file = open(os.path.join(this_directory, 'ellis_island.fasta'), 'w')
        # first get the unique members in these families, and make a list of their accessions
        family_list = [fam.family.id for fam in self.get_families_passing_stage1()]
        uniprot_set = find_uniprots_for_family_list(family_list)
        # print str(len(uniprot_set)) + ' uniprot ids in stage 1 families'

        # now we need to compile the FASTA for these uniprot sequences
        queue = Queue.Queue()
        for i in range(THREADS):
            t = threading.Thread(target=self.compile_fasta_for_ellis_island, args=(queue, ))
            t.setDaemon(True)
            t.start()
        for chunk in chunkify(uniprot_set, THREADS):
            if chunk:
                queue.put((chunk, this_directory), True)
        queue.join() # wait for completion
        temporary_fasta_paths = [file for file in glob.glob(os.path.join(this_directory,'job%d*'%self.id))]
        # concatenate these member fasta files
        args = ['cat'] + temporary_fasta_paths
        process = subprocess.Popen(args, shell=False, stdout=fasta_file, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        process.communicate()
        fasta_file_path = fasta_file.name
        fasta_file.close()
        
        # Now, get the top k1 sequences to the query
        # first, get an hmm of the query
        hmm_file = open(tempfile.mkstemp(dir=this_directory)[1],'w')
        hmm_file.write(hmmbuild(self.fasta, 'query')[0])
        hmm_file.close()
        ((full_length, domain), error) = hmmsearch(hmm_file.name, fasta_file_path, 
                                    options=['--cpu', '8', '-E', '13-3', '--domE','1e-3'])

        # get a list of the accessions sorted by bit score and only those that pass coverage conditions
        # this already has an implicit evalue threshold imposed by hmmer
        # DD: added reverse = true as was sorted ascending before
        top_scores = [(row.target_name.split('|')[1], row.full_score) for row in sorted(parse_domtblout(domain), 
            key=lambda hmmerrow: hmmerrow.full_score, reverse=True) 
            if float(row.ali_to - row.ali_from)/len(self.fasta_sequence) >= self.ellis_island_member_coverage and
            float(row.hmm_to - row.hmm_from)/int(row.target_length) >= self.ellis_island_query_coverage]
        # print str(len(top_scores)) + ' seqs before taxon clustering'

        # group by taxon id 
        top_grouped_by_taxon = []
        for entry in top_scores:
            u = UniProt.objects.get(accession = entry[0])
            top_grouped_by_taxon.append([u.taxon, u, entry[1]]);
        list.sort(top_grouped_by_taxon, key=lambda hit: hit[0].id) # must sort before groupby
        
        # for each taxon
        passing_accessions = []
        for taxon, taxon_members in groupby(top_grouped_by_taxon, key=lambda hit: hit[0]):
            taxon_members_list = list(taxon_members) # gonna use twice so iterator won't do
            if len(taxon_members_list) == 1: # only 1 seq - no need to cluster
                passing_accessions.append((taxon_members_list[0][1].accession, taxon_members_list[0][2]))
            else: # call uclust
                taxon_fasta = ''
                best_bitscore = 0.0
                members_dict = {}
                member_count = 0
                for m in taxon_members_list:
                    member_count += 1
                    taxon_fasta += m[1].fasta_from_db + '\n'
                    members_dict[m[1].accession] = m[2]
                    if m[2] > best_bitscore:
                        best_bitscore = m[2]
                        highest_scoring_member = m
                cluster_membership = uclust(taxon_fasta, this_directory, self.ellis_island_taxon_clustering_pid)

                # get cluster with highest (bit)scoring member
                highest_scoring_cluster = -1
                for member in cluster_membership:
                    if member[1] == highest_scoring_member[1].accession:
                        highest_scoring_cluster = member[0]
                for member in cluster_membership:
                    if member[0] == highest_scoring_cluster:
                        passing_accessions.append((member[1], members_dict[member[1]]))
        # print str(len(passing_accessions)) + ' seqs after taxon clustering'

        # sort by bit score, get the best, then chop off the bitscore column
        list.sort(passing_accessions, key=lambda x: x[1], reverse=True)
        top_accessions = [item[0] for item in passing_accessions[0:self.ellis_island_maxseqs]] 
        
        family_list = find_families_for_uniprot_accession_list(top_accessions, family_list)
        families_needed = []
        taxa_covered = 0
        while top_accessions:
            this_family = ellis_island_get_most_taxa(top_accessions, [fam.id for fam in family_list])
            for member in get_members_for_elimination(top_accessions, this_family):
                top_accessions.remove(member)
                
            family_list.remove(this_family)
            families_needed.append(this_family)
        
        # change permissions of everything in this directory, and remove it
        # not working yet (leaves empty directory)
        subprocess.call(['chmod', '-R', '777', this_directory])
        subprocess.call(['rm', '-rf', this_directory])
        #shutil.rmtree(this_directory)
        #subprocess.call(['rmdir', this_directory])

        return families_needed[0:self.ellis_island_maxfams]

    # thread function called by scan_subtree_hmms_mt
    def threaded_subtree_hmm(self, queue, exec_lock):
        while True:
            fam = queue.get()
            #print str(fam.family) + ' ' + str(len(fam.family.canonical_root_node().get_included_leaves()))
            fam.subtree_hmm_scan_mt(exec_lock)
            #fam.save()
            connection.close()
            queue.task_done()
    
    # multithreaded version of scan_subtree_hmms()
    def scan_subtree_hmms_mt(self, stage_2_families):
        queue = Queue.Queue()
        exec_lock = threading.Lock() # hmmscan instances have to start one at a time
        for i in range(THREADS):
            t = threading.Thread(target=self.threaded_subtree_hmm, args=(queue, exec_lock, ))
            t.setDaemon(True)
            t.start()
            
        for fam in stage_2_families:
            if fam.passes_stage1_coverage_conditions:
                queue.put(fam, True)
            else:
                # This is kind of stupid.  Should never really be used 
                # for families that don't pass stage 1 coverage.  
                # They are never scanned in stage 2.
                fam.successful_hmmscan = False
                fam.save()
        queue.join()
        
    def scan_subtree_hmms(self, stage_2_families):
        for fam in stage_2_families:
            if fam.passes_stage1_coverage_conditions:
                if fam.subtree_hmm_scan():
                    fam.successful_hmmscan = True
                    fam.save()
                else:
                    fam.successful_hmmscan = False
                    fam.save()
            else:
                # This is kind of stupid.  Should never really be used 
                # for families that don't pass stage 1 coverage.  
                # They are never scanned in stage 2.
                fam.successful_hmmscan = False
                fam.save()
        return None

    def analyze_top_scoring_nodes(self):
        for fam in self.families.all():
            if fam.passes_stage1_coverage_conditions:
                fam.passed_stage1 = True
                fam.save()
                if fam.successful_hmmscan:
                    if fam.passes_stage2_coverage_conditions:
                        if fam.get_enclosing_clade_root() is not None:
                            fam.enclosing_clade_root = fam.get_enclosing_clade_root()
                            ali = fam.top_scoring_node.alignment_to_query(self.fasta)
                            fam.query_hmm_consensus_pid = ali[0]
                            fam.query_hmm_alignment = repr(ali[1])                            
                            # get description of enclosing clade
                            try:
                                fam.enclosing_clade_description = fam.enclosing_clade_root.consensus_description.get().consensus_uniprot_description
                            except:
                                fam.enclosing_clade_description = fam.enclosing_clade_root.get_treenode_names()['treenode_name']
                            # get mrca of enclosing clade
                            try:
                                fam.enclosing_clade_taxonomy = fam.enclosing_clade_root.mrca.get().mrca
                            except:
                                fam.enclosing_clade_taxonomy = fam.enclosing_clade_root.get_taxon()
                            fam.save()
                            # does the family pass stage 2 query to hmm criteria?
                            if fam.query_hmm_consensus_pid >= self.stage2_query_hmm_pid:
                                # if yes, this family passed stage 2
                                fam.passed_stage2 = True
                                fam.reason_eliminated = 'None'
                                fam.save()
                            else:
                                fam.reason_eliminated = 'Insufficient similarity between query and HMM consensus of top-scoring node'
                                fam.save()
                        else:
                            fam.reason_eliminated = 'No Enclosing Clade'
                            fam.save()
                    else:
                        fam.reason_eliminated = 'Insufficient coverage between query and HMM consensus of top-scoring node'
                        fam.save()
                else:
                    fam.reason_eliminated = 'hmmscan error in stage 2'
                    fam.save()
            else:
                fam.reason_eliminated = 'Insufficient coverage between query and HMM consensus of family root node'
                fam.save()
        return None
    
    def get_member_fasta(self):
        # gets the fasta for the members
        members = set()
        fasta = ''
        seq_count = 0
        seqs_found = False
        for family in self.get_families_passing_stage2():
            for leaf in family.enclosing_clade_root.get_included_leaves():
                if ((leaf.sequence_header) and (leaf.sequence_header.uniprot) and
                    (leaf.sequence_header.uniprot.fasta_from_db) and 
                    (leaf.sequence_header.uniprot not in members)):
                    seq_count += 1
                    seqs_found = True
                    members.add(leaf.sequence_header.uniprot)
                    fasta += leaf.sequence_header.uniprot.fasta_from_db + '\n'
        #print 'member fasta size pre k-tuple: ' + str(seq_count)
        return fasta 

    def build_fatcat_job_msa(self):
        # builds an msa from the k closest sequences to the query
        # returns a biopython alignio object
        # first get the query fasta
        member_fasta = self.get_member_fasta()
        if not member_fasta:
            return None
            
        member_query_fasta = '>query' + self.fasta_header + '\n' + self.fasta_sequence + '\n' + member_fasta

        # now get the k closest and construct the alignment
        closest = get_top_k_to_query(member_query_fasta, self.stage2_max_sequences)
        
        # align these with mafft and return
        return AlignIO.read(StringIO(mafft(closest)), 'fasta')
        # return AlignIO.read(StringIO(mafft(closest, ['--retree', '1'])), 'fasta') # test

    def fetch_alignment(self):
        alignment = ''
        for record in self.members.filter(is_query = True):
            alignment += '>QUERY\n'
            alignment += record.alignment + '\n'
        for record in self.members.filter(is_query = False):
            alignment += '>' + record.uniprot.uniprot_identifier + ' ' + record.uniprot.description + '\n'
            alignment += record.alignment + '\n'
        return alignment
            
    def build_fatcat_job_tree(self, msa):
        # builds the tree for this fatcat job.  takes the alignio object from above
        # and masks the tree accordingly.  returns the newick string
        masked_msa = chop_msa(msa, 0.7)
        # build the full tree of k sequences
        return fast_tree_mp(masked_msa)
        
    def cut_fatcat_job_tree(self, msa, tree):
        # do a tree cut on this and return...we should maybe just redo this part?
        # msa is an alignio object, tree is a newick string, so it has to be read in
        # this is totally stupid, we need to rewrite these things so we're not going
        # from one piece of thing to another so often.
        tree_cut = MSATree(msa, Phylo.read(StringIO(tree),'newick'))
        [sub_tree, sub_msa] = tree_cut.maximal_subtree_msa('query' + self.fasta_header.split()[0], self.stage3_max_sequences)
        # save the newick string
        self.tree_string = sub_tree.format('newick')
        self.phyloxml_string = sub_tree.format('phyloxml')
        self.save()
        return [sub_tree, sub_msa]
    


    def store_functions_transferred(self, gene=None, description=None, ec=None):
        json_object = {}
        if gene:
            json_object['gene'] = gene
        if description:
            json_object['description'] = description
        if ec:
            json_object['ec'] = ec
        self.functions = json.dumps(json_object)
        self.save()
        return None

    def get_pathologic_string(self, gene, description, ec, unique):
        if not gene:
            return None
        # todo: more elegant way on ensuring no duplicates
        pathologic = 'ID\t' + gene[0][0] + '_' + str(unique) + '\n'
        pathologic += 'NAME\t' + gene[0][0] + '\n'
        
        # take the top function
        # and take functions after that if their weight is greater than half the one before it
        # idea is to give more than one if the weights are closely matched
        # this could be replaced at some point
        if description:
            prev_weight = 0
            for d in description:
                if d[1] > prev_weight / 2:
                    pathologic += 'FUNCTION\t' + d[0] + '\n'
                else:
                    break
                prev_weight = d[1]
    
        if ec:
            pathologic += 'EC\t' + ec[0][0] + '\n'
        
        # type = protein plus double slash to end
        pathologic += 'PRODUCT-TYPE\tP\n//\n'
        
        return pathologic
    
    def predict_query_ec_number(self, base=0.5):
        ec_names = {}
        ec_list = []
        total_weighted_votes = 0
        for ortholog in self.get_orthologs():
            query_pid = ortholog.member_to_query_pid
            weight = pow(base, (10-(query_pid//10)))
            try:
                this_ec = UniProtEC.objects.get(uniprot=ortholog.uniprot)
                this_ec_string = str(this_ec.ec)
            except:
                continue
            if this_ec.is_in_brenda_f:
                # Djesus, why 5??!?
                ec_names[this_ec_string] = ec_names.get(this_ec_string, 0) + 5*weight
                total_weighted_votes += 5*weight
            else:
                ec_names[this_ec_string] = ec_names.get(this_ec_string, 0) + weight
                total_weighted_votes += weight
        for ec, votes in ec_names.items():
            ec_list.append( ( ec, float(1.0*votes/total_weighted_votes) ) )     
        return sorted(ec_list, key=lambda g: g[1], reverse=True)
    
    def predict_query_gene_name(self, base=0.5):
        gene_names = {}
        gene_list = []
        total_weighted_votes = 0
        for ortholog in self.get_orthologs():
            query_pid = ortholog.member_to_query_pid
            weight = pow(base, (10-(query_pid//10)))
            genes = UniProtGene.objects.filter(uniprot = ortholog.uniprot)
            if genes:
                gene = genes[0].name.upper()
                if ortholog.uniprot.in_swissprot_f:
                    gene_names[gene] = gene_names.get(gene, 0) + 10*weight
                    total_weighted_votes += 10*weight
                else:
                    gene_names[gene] = gene_names.get(gene, 0) + weight
                    total_weighted_votes += weight
        for gene, votes in gene_names.items():
            gene_list.append( ( gene, float(1.0*votes/total_weighted_votes) ) )     
        return sorted(gene_list, key=lambda g: g[1], reverse=True)

    def predict_query_description(self, base=0.5):
        descriptions = {}
        description_list = []
        total_weighted_votes = 0
        for ortholog in self.get_orthologs():
            description = ortholog.uniprot.description
            query_pid = ortholog.member_to_query_pid
            weight = pow(base, (10-(query_pid//10)))
            if ortholog.uniprot.in_swissprot_f:
                descriptions[description] = descriptions.get(description, 0) + 10*weight
                total_weighted_votes += 10*weight
            else:
                descriptions[description] = descriptions.get(description, 0) + weight
                total_weighted_votes += weight
        for d, v in descriptions.items():
            description_list.append( ( d, float(1.0*v/total_weighted_votes) ) )
        return sorted(description_list, key=lambda desc: desc[1], reverse=True)

    def get_fatcat_job_msa(self, masking=None):
        msa = ''
        for member in self.members.all():
            record = member.get_full_alignment()
            if not record.endswith('\n'):
                record += '\n'
            msa += record
        return msa

    def send_email(self):
        # sends email
        if self.user_email:
            text_body = '''PhyloFacts FAT-CAT Job %d has completed.
                \n\nhttp://makana.berkeley.edu/fatcat/%d/
                \n%s
                \n\n%s\n%s
                \n\n''' % (self.id, self.id, self.fatcat_job_name, self.fasta_header, self.fasta_sequence[:10])
            html_body = '''<h2>PhyloFacts FAT-CAT Job %d has completed</h2>.
                \n\n<p><a href="http://makana.berkeley.edu/fatcat/%s">click here to view the results</a></p>
                \n\n<p>%s</p>
                \n\n<p>%s</p>\n<p>%s</p>
                \n\n''' % (self.id, self.id, self.fatcat_job_name, self.fasta_header, self.fasta_sequence[:10])
            PhyloFactsEmail.objects.create(
                recipient = self.user_email,
                text_body = text_body,
                html_body = html_body,
                subject = self.email_subject
            )
        return None        

    def get_orthologs(self):
        # this function should return a list of the orthologs for this job.
        # this will probably be changed many, many times
        return list(self.members.filter(classification='Orthologs'))

    def get_ortholog_file(self, delimiter='\t'):
        # builds the orthologs downloadable file for this job.
        fields = ['Uniprot Accession', 'Species', 'Taxon ID', 'Description', '% ID to query', 'Query coverage', 'Hit coverage']
        file = delimiter.join(fields) + '\n' 
        for ortholog in self.get_orthologs():
            file += '"%s"' % (ortholog.uniprot.accession)
            file += delimiter + '"%s"' % (ortholog.uniprot.taxon.scientific_name)
            file += delimiter + '%d' % (ortholog.uniprot.taxon.id)
            file += delimiter + '"%s"' % (ortholog.uniprot.description)
            file += delimiter + '%.2f' % (ortholog.member_to_query_pid)
            file += delimiter + '%.2f' % (ortholog.query_coverage)
            file += delimiter + '%.2f\n' % (ortholog.member_coverage)
        return file

    def populate_member_go_annotations(self):
        ''' This function populates the members with their go annotations '''
        for member in self.get_orthologs():
            annotations = {}
            for annotation in  Annotation.objects.filter(uniprot__accession = member.uniprot.accession, 
                    ontology_term__ontology = 'go').values('id','ontology_term__accession', 'ontology_term__name', 'evidence__code', 'evidence__name', 'evidence__priority'):
                if annotation['ontology_term__accession'] not in annotations:
                    annotations[annotation['ontology_term__accession']] = annotation
                elif annotation['evidence__priority'] < annotations[annotation['ontology_term__accession']]['evidence__priority']:
                    annotations[annotation['ontology_term__accession']] = annotation
            # populate the go annotations
            for (accession, annotation) in annotations.items():
                member.go_annotations.create(annotation_id = annotation['id'], fatcat_job_id=self.id)
        return None

    def get_summary_go_annotations(self):
        ''' This will retrieve the already populated go annotations.  The dict returned will have three keys
            one for molecular function, one for biological process, and one for cellular component. 
            these dicts will contains lists of go annotations that are found in the job's orthologs. '''
        try:
            return self._get_summary_go_annotations
        except:
            pass

        go_annotations = {}
        go_annotations['molecular_function'] = {}
        go_annotations['biological_process'] = {}
        go_annotations['cellular_component'] = {}
        for member in self.get_orthologs():
            these_annotations = member.go_annotations.all()
            if these_annotations:
                for this_annotation in these_annotations:
                    subset = OntologySubsetMembership.objects.filter(term=this_annotation.annotation.ontology_term)[0].subset.name
                    this_go_annotation = {
                        'accession': this_annotation.annotation.ontology_term.accession,
                        'description':  this_annotation.annotation.ontology_term.name,
                        'evidence_code':  this_annotation.annotation.evidence.code,
                        'evidence_name':  this_annotation.annotation.evidence.name,
                        'evidence_priority':  this_annotation.annotation.evidence.priority
                    }
                    if this_go_annotation['accession'] in go_annotations[subset]: 
                        if (this_go_annotation['evidence_priority'] < go_annotations[subset][this_go_annotation['accession']]['evidence_priority']):                 
                            go_annotations[subset][this_go_annotation['accession']] = this_go_annotation
                    else:
                        go_annotations[subset][this_go_annotation['accession']] = this_go_annotation
        return_dict = {}
        return_dict['molecular_function'] = sorted(go_annotations['molecular_function'].values(), key=lambda annotation: annotation['evidence_priority'])
        return_dict['cellular_component'] = sorted(go_annotations['cellular_component'].values(), key=lambda annotation: annotation['evidence_priority'])
        return_dict['biological_process'] = sorted(go_annotations['biological_process'].values(), key=lambda annotation: annotation['evidence_priority'])
        self._get_summary_go_annotations = return_dict
        return self._get_summary_go_annotations

    def log_job_completion(self):
        now = datetime.now()
        td = self.created_at - now
        self.job_completion_time = float(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
        self.is_done = True
        self.save()
        self.send_email()
        return None
        
    class Meta:
        db_table = u'fatcat2_job'

class FatcatJobFamily(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name='families')
    family = models.ForeignKey(Family, related_name='fatcat_pf_family')
    family_description = models.TextField()
    family_taxonomy = models.ForeignKey(UniProtTaxonomy, related_name='fam_taxonomy')
    family_mda = models.TextField()
    family_type = models.TextField()
    family_e_value = models.FloatField()
    top_scoring_node = models.ForeignKey(TreeNode, related_name='tsn')
    top_scoring_node_e_value = models.FloatField()
    top_scoring_node_description = models.TextField()
    top_scoring_node_taxonomy = models.ForeignKey(UniProtTaxonomy, related_name='tsn_taxonomy')
    top_scoring_node_mda = models.TextField()
    enclosing_clade_root = models.ForeignKey(TreeNode, related_name='ec')
    enclosing_clade_description = models.TextField()
    enclosing_clade_taxonomy = models.ForeignKey(UniProtTaxonomy, related_name='ec_taxonomy')
    enclosing_clade_mda = models.TextField()
    query_hmm_consensus_pid = models.FloatField()
    query_hmm_alignment = models.TextField()
    successful_hmmscan = models.BooleanField()
    reason_eliminated = models.TextField()
    passed_stage1 = models.BooleanField()
    passed_stage2 = models.BooleanField()

    @property
    def passes_stage1_coverage_conditions(self):
        md = self.match_data.get()
        if self.family_type == 'Pfam':
            return (md.stage1_hmm_coverage >= self.fatcat_job.stage1_pfam_family_hmm_coverage)
        else:
            return ((md.stage1_hmm_coverage >= self.fatcat_job.stage1_mda_family_hmm_coverage) and
                    (md.stage1_query_coverage >= self.fatcat_job.stage1_mda_family_query_coverage))

    @property
    def passes_stage2_coverage_conditions(self):
        md = self.match_data.get()
        if self.family_type == 'Pfam':
            return (md.stage2_hmm_coverage >= self.fatcat_job.stage2_pfam_family_hmm_coverage)
        else:
            return ((md.stage2_hmm_coverage >= self.fatcat_job.stage2_mda_family_hmm_coverage) and
                    (md.stage2_query_coverage >= self.fatcat_job.stage2_mda_family_query_coverage))

    @property
    def has_enclosing_clade(self):
        return self.enclosing_clade_root is not None

    @property
    def member_fasta_path(self):
        return '/clusterfs/ohana/software/fatcat/phylofacts_member_fasta/%s.fasta' % self.family.get_accession()

    def get_family_link(self):
        acc = self.family.get_accession()
        return '<a target="_blank" href="/phylofacts/family/%s/">%s</a>' % (acc, acc)

    def get_family_taxonomy_link(self):
        if (self.family_taxonomy is not None):
            tax = self.family_taxonomy
            if (tax.id < 2):
                if (tax.common_name):
                    return '<span class="tip" title="%s">%s</span>' % (tax.common_name, tax.scientific_name)
                else:
                    return '<span>%s</span>' % (tax.scientific_name)
            else:
                if (tax.common_name):
                    return '<a class="tip" target="_blank" title="%s" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.common_name, tax.id, tax.scientific_name)
                else:
                    return '<a target="_blank" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.id, tax.scientific_name)              
        else:
            return ''

    def get_enclosing_clade_taxonomy_link(self):
        if (self.enclosing_clade_taxonomy is not None):
            tax = self.enclosing_clade_taxonomy
            if (tax.id < 2):
                if (tax.common_name):
                    return '<span class="tip" title="%s">%s</span>' % (tax.common_name, tax.scientific_name)
                else:
                    return '<span>%s</span>' % (tax.scientific_name)
            else:
                if (tax.common_name):
                    return '<a class="tip" target="_blank" title="%s" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.common_name, tax.id, tax.scientific_name)
                else:
                    return '<a target="_blank" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.id, tax.scientific_name)              
        else:
            return ''
   
    def get_top_scoring_node_taxonomy_link(self):
        if (self.top_scoring_node_taxonomy is not None):
            tax = self.top_scoring_node_taxonomy
            if (tax.id < 2):
                if (tax.common_name):
                    return '<span class="tip" title="%s">%s</span>' % (tax.common_name, tax.scientific_name)
                else:
                    return '<span>%s</span>' % (tax.scientific_name)
            else:
                if (tax.common_name):
                    return '<a class="tip" target="_blank" title="%s" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.common_name, tax.id, tax.scientific_name)
                else:
                    return '<a target="_blank" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.id, tax.scientific_name)              
        else:
            return ''
    
    def get_distant_clade_taxonomy_link(self):
        if (self.top_scoring_node_taxonomy is not None):
            tax = self.top_scoring_node_taxonomy
            if (tax.id < 2):
                if (tax.common_name):
                    return '<span class="tip" title="%s">%s</span>' % (tax.common_name, tax.scientific_name)
                else:
                    return '<span>%s</span>' % (tax.scientific_name)
            else:
                if (tax.common_name):
                    return '<a class="tip" target="_blank" title="%s" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.common_name, tax.id, tax.scientific_name)
                else:
                    return '<a target="_blank" href="http://www.uniprot.org/taxonomy/%d">%s</a>' % (tax.id, tax.scientific_name)              
        else:
            return ''

    def get_query_hmm_pid_link(self): 
        if self.query_hmm_consensus_pid is not None:
            ret_string = '<a href="javascript:;" onclick="startJalview(' 
            ret_string += "%s, 'subtree alignment', 'alwvr_%d')" % (self.query_hmm_alignment, self.family.id) 
            ret_string += '">%.1f</a>' % (self.query_hmm_consensus_pid)
            return ret_string
        else:
            return '<span style="color:red">N/A</span>'

    def get_enclosing_clade_support_links(self):
        ret_string = '<div class="btn-toolbar ec-support-links"><div class="btn-group">'
        ret_string += '<a class="btn" href="javascript:;" onclick="window.open('
        ret_string += "'/fatcat/family/%d/tree/', 'height=800,width=800,left=100,top=100,resizable=yes,scrollbars=yes,toolbar=yes,menubar=no,location=no,directories=no, status=yes');" % (self.id)
        ret_string += '"><img class="tree-image tip" title="View the Enclosing Clade" src="/static/img/icons/tree-icon.png"></a><a class="btn tip" href="/phylofacts/tree_node_view/%d/" title="View more information about the top-scoring node" target="_blank"><img class="annotation-image" src="/static/img/icons/annotation-icon.png"></a></div></div>' % (self.top_scoring_node.id)
        return ret_string

    def get_distant_clade_support_links(self):
        ret_string = '<div class="btn-toolbar dc-support-links"><div class="btn-group">'
        ret_string += '<a class="btn" href="javascript:;" onclick="window.open('
        ret_string += "'/phylofacts/tree_node_view/%d/orthoscope/', 'height=800,width=800,left=100,top=100,resizable=yes,scrollbars=yes,toolbar=yes,menubar=no,location=no,directories=no, status=yes');" % (self.top_scoring_node.id)
        ret_string += '"><img class="tree-image tip" title="View the Enclosing Clade" src="/static/img/icons/tree-icon.png"></a><a class="btn tip" href="/phylofacts/tree_node_view/%d/" target="_blank" title="View more information about the top-scoring node"><img class="annotation-image"  src="/static/img/icons/annotation-icon.png"></a></div></div>' % (self.top_scoring_node.id)
        return ret_string

    def get_reason_eliminated_link(self):
        return '<div style="color:red" class="tip" title="%s">?</div>' % self.reason_eliminated

    def get_family_row_datatable_object(self):
        try:
            return self._get_family_row_datatable_object
        except:
            pass
        
        if self.family_type == 'Pfam':
            type = 'Pfam'
        else:
            type = 'MDA'
        md = self.match_data.get()
        self._get_family_row_datatable_object = {
            'DT_RowId': 'fam_%d' % self.family.id,
            '0': self.get_family_link(),
            '1': type,
            '2': self.family_description,
            '3': self.get_family_taxonomy_link(),
            '4': self.family_e_value,
            '5': '%d - %d' % (md.stage1_ali_from, md.stage1_ali_to),
            '6': '%.1f' % (100 * md.stage1_hmm_coverage),
            '7': '%.1f' % (100 * md.stage1_query_coverage)
        }
        return self._get_family_row_datatable_object
             
    def get_enclosing_clade_row_datatable_object(self):
        try:
            return self._get_enclosing_clade_row_datatable_object
        except:
            pass
        
        if self.family_type == 'Pfam':
            type = 'Pfam'
        else:
            type = 'MDA'
        md = self.match_data.get()
        self._get_enclosing_clade_row_datatable_object = {
            'DT_RowId': 'ec_%d' % self.family.id,
            '0': self.get_family_link(),
            '1': type,
            '2': self.enclosing_clade_description,
            '3': self.get_top_scoring_node_taxonomy_link(),
            '4': self.get_query_hmm_pid_link(),
            '5': '%.1f' % (100 * md.stage2_query_coverage),
            '6': '%.1f' % (100 * md.stage2_hmm_coverage),
            '7': '%d - %d' % (md.stage2_ali_from, md.stage2_ali_to),
            '8': self.top_scoring_node_e_value,
            '9': self.get_enclosing_clade_support_links()            
        }
        return self._get_enclosing_clade_row_datatable_object
    
    def get_distant_clade_row_datatable_object(self):
        try:
            return self._get_distant_clade_row_datatable_object
        except:
            pass
        
        if self.family_type == 'Pfam':
            type = 'Pfam'
        else:
            type = 'MDA'
        md = self.match_data.get()
        self._get_distant_clade_row_datatable_object = {
            'DT_RowId': 'dc_%d' % self.family.id,
            '0': self.get_reason_eliminated_link(),
            '1': self.get_family_link(),
            '2': type,
            '3': self.top_scoring_node_description,
            '4': self.get_distant_clade_taxonomy_link(),
            '5': self.top_scoring_node_e_value,
            '6': self.get_query_hmm_pid_link(),
            '7': '%.1f' % (100 * md.stage2_query_coverage),
            '8': '%.1f' % (100 * md.stage2_hmm_coverage),
            '9': '%d - %d' % (md.stage2_ali_from, md.stage2_ali_to),
            '10': self.get_distant_clade_support_links()            
        }
        return self._get_distant_clade_row_datatable_object
    
    def get_enclosing_clade_root(self):
        if (not self.passes_stage1_coverage_conditions) or (not self.successful_hmmscan):
            return None

        # get enclosing clade root here
        orthology_methods = self.fatcat_job.stage2_orthology_methods.strip().split()
        orthology_methods = [x + '_supported_ancestor' for x in orthology_methods]
        num_orthology_methods_required = self.fatcat_job.stage2_redundant_methods_required
        kt = self.fatcat_job.stage2_enclosing_clade_kerf_threshold
        orthology_items = [x for x in self.top_scoring_node.orthology_supported_ancestors(kerf_threshold=kt).items() 
                            if x[1] is not None and x[0] in orthology_methods]
        if orthology_items:
            # sort by left id, nodes higher up the tree appear first
            orthology_items = sorted(orthology_items, key=lambda t: t[1].left_id)
        else:
            return None
        matched = 0
        top_enclosing_clade = None
        for (system, root) in orthology_items:
            if root == top_enclosing_clade:
                matched += 1
            else:
                top_enclosing_clade = root
                matched = 1
            if matched == num_orthology_methods_required:
                return top_enclosing_clade
        return None 
    
    def subtree_hmm_scan(self):
        # try-excepts like these are hard to debug
        try:
            #hmm_file = self.family.subtree_hmm_file(no_build = True)
            hmm_file = self.family.subtree_hmm_file()
            if hmm_file is None:
                return False
            # hmmscan options for different parameters
            # we only give one cpu as we'll be running lots of 
            if self.fatcat_job.parameterization == 'High recall':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'High precision':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'Remote homologs':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'Partial sequence':
                hmmscan_options = ['--max', '-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'Custom':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            else:
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]

            ((normal, domain), error_code) = hmmscan(self.fatcat_job.fasta, hmm_file, hmmscan_options)
            if error_code:
                return False
            domain_table = parse_domtblout(domain)
            domain_table.sort(reverse=True)
            if (len(domain_table) > 0):
                tsn = TreeNode.objects.get(id=int(domain_table[0].target_name))
                self.top_scoring_node = tsn
                self.top_scoring_node_e_value = domain_table[0].this_domain_i_value
                # get top scoring node description
                try:
                    self.top_scoring_node_description = tsn.consensus_description.get().consensus_uniprot_description
                except:
                    self.top_scoring_node_description = tsn.get_treenode_names()['treenode_name']
                # get top scoring node taxonomy
                try:
                    self.top_scoring_node_taxonomy = tsn.mrca.get().mrca
                except:
                    self.top_scoring_node_taxonomy = tsn.get_taxon()
                # get top scoring node consensus mda
                '''try:
                    self.top_scoring_node_mda = 
                except:
                    self.top_scoring_node_mda = '''
                self.save()
                md = self.match_data.get()
                md.stage2_hmm_from = domain_table[0].hmm_from
                md.stage2_hmm_to = domain_table[0].hmm_to
                md.stage2_ali_from = domain_table[0].ali_from
                md.stage2_ali_to = domain_table[0].ali_to
                md.stage2_hmm_length = domain_table[0].target_length
                md.stage2_region_aligned = md.stage2_ali_to - md.stage2_ali_from
                md.stage2_query_coverage = float(md.stage2_ali_to - md.stage2_ali_from)/len(self.fatcat_job.fasta_sequence)
                md.stage2_hmm_coverage = float(md.stage2_hmm_to - md.stage2_hmm_from)/int(md.stage2_hmm_length)
                md.save()
                return True
            else:
                return False
        except:
            return False

    def hmmscan_mt(self, exec_lock, fasta, hmm_file, options=['-E', '1E-3']):
        #TODO input validation on E
        #args = ["/clusterfs/ohana/home/ddineen/stuff/hmmer-3.0/src/hmmscan", 
        args = ["/clusterfs/ohana/software/bin/hmmscan", 
            "--cpu", "1",
            "-o", "/dev/null", 
            "--notextw", 
            "--tblout", "/dev/stdout", 
            "--domtblout", "/dev/stderr",
            ] + options + [ 
            hmm_file,
            "/dev/stdin" ]

        exec_lock.acquire()
        process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        exec_lock.release()
        return (process.communicate(fasta), process.returncode)

    def subtree_hmm_scan_mt(self, exec_lock):
        # try-excepts like these are hard to debug
        try:
            #hmm_file = self.family.subtree_hmm_file(no_build = True)
            hmm_file = self.family.subtree_hmm_file()
            if hmm_file is None:
                return False
            # hmmscan options for different parameters
            # we only give one cpu as we'll be running lots of 
            if self.fatcat_job.parameterization == 'High recall':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'High precision':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'Remote homologs':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'Partial sequence':
                hmmscan_options = ['--max', '-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            elif self.fatcat_job.parameterization == 'Custom':
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]
            else:
                hmmscan_options = ['-E', str(self.fatcat_job.stage2_e_value), '--domE', str(self.fatcat_job.stage2_e_value)]

            #((normal, domain), error_code) = hmmscan(self.fatcat_job.fasta, hmm_file, hmmscan_options)
            ((normal, domain), error_code) = self.hmmscan_mt(exec_lock, self.fatcat_job.fasta, hmm_file, hmmscan_options)
            if error_code:
                self.successful_hmmscan = False
                self.save()
                return False
            domain_table = parse_domtblout(domain)
            domain_table.sort(reverse=True)
            if (len(domain_table) > 0):
                tsn = TreeNode.objects.get(id=int(domain_table[0].target_name))
                self.top_scoring_node = tsn
                self.top_scoring_node_e_value = domain_table[0].this_domain_i_value
                # get top scoring node description
                try:
                    self.top_scoring_node_description = tsn.consensus_description.get().consensus_uniprot_description
                except:
                    self.top_scoring_node_description = tsn.get_treenode_names()['treenode_name']
                # get top scoring node taxonomy
                try:
                    self.top_scoring_node_taxonomy = tsn.mrca.get().mrca
                except:
                    self.top_scoring_node_taxonomy = tsn.get_taxon()
                # get top scoring node consensus mda
                '''try:
                    self.top_scoring_node_mda = 
                except:
                    self.top_scoring_node_mda = '''
                self.save()
                md = self.match_data.get()
                md.stage2_hmm_from = domain_table[0].hmm_from
                md.stage2_hmm_to = domain_table[0].hmm_to
                md.stage2_ali_from = domain_table[0].ali_from
                md.stage2_ali_to = domain_table[0].ali_to
                md.stage2_hmm_length = domain_table[0].target_length
                md.stage2_region_aligned = md.stage2_ali_to - md.stage2_ali_from
                md.stage2_query_coverage = float(md.stage2_ali_to - md.stage2_ali_from)/len(self.fatcat_job.fasta_sequence)
                md.stage2_hmm_coverage = float(md.stage2_hmm_to - md.stage2_hmm_from)/int(md.stage2_hmm_length)
                md.save()
                self.successful_hmmscan = True
                self.save()
                return True
            else:
                self.successful_hmmscan = False
                self.save()
                return False
        except:
        #    print 'failed'
        #    print sys.last_traceback
            self.successful_hmmscan = False
            self.save()
            return False

    class Meta:
        db_table = u'fatcat2_job_family'

class FatcatJobFamilyMatchData(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job_family = models.ForeignKey(FatcatJobFamily, related_name='match_data')
    stage1_ali_from = models.IntegerField()
    stage1_ali_to = models.IntegerField()
    stage1_hmm_from = models.IntegerField()
    stage1_hmm_to = models.IntegerField()
    stage1_hmm_length = models.IntegerField()
    stage1_hmm_coverage = models.FloatField()
    stage1_query_coverage = models.FloatField()
    stage1_region_aligned = models.IntegerField()
    stage2_hmm_from = models.IntegerField()
    stage2_hmm_to = models.IntegerField()
    stage2_ali_from = models.IntegerField()
    stage2_ali_to = models.IntegerField()
    stage2_hmm_coverage = models.FloatField()
    stage2_query_coverage = models.FloatField()
    stage2_region_aligned = models.IntegerField()
    stage2_hmm_length = models.IntegerField()

    class Meta:
        db_table = u'fatcat2_job_family_match_data'    

class FatcatJobPfams(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name="pfams")
    pfam_description = models.TextField()
    pfam_accession = models.TextField()
    pfam_shortname = models.TextField()
    ali_from = models.IntegerField()
    ali_to = models.IntegerField()
    hmm_from = models.IntegerField()
    hmm_to = models.IntegerField()
    c_evalue = models.FloatField()
    i_evalue = models.FloatField()

    class Meta:
        db_table = u'fatcat2_job_pfams'

class FatcatJobMemberClusters(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name='clusters')
    cluster_representative = models.ForeignKey('FatcatJobMembers')
    taxon = models.ForeignKey(UniProtTaxonomy)
    cluster_size = models.IntegerField()
    
    def get_uniprot_id_link(self):
        return '<a href="http://www.uniprot.org/uniprot/%s" target="_blank">%s</a>' % (self.cluster_representative.uniprot.accession, self.cluster_representative.uniprot.uniprot_identifier)

    def get_cluster_taxonomy_link(self):
        if self.taxon.common_name:
            ret_string = '<table class="cluster-taxon-table"><tr><td><a href="http://www.uniprot.org/taxonomy/%d" class="tip" title="%s" target="_blank">%s</a></td>' % (self.taxon.id, self.taxon.common_name, self.taxon.scientific_name)
        else:
            ret_string = '<table class="cluster-taxon-table"><tr><td><a href="http://www.uniprot.org/taxonomy/%d" target="_blank">%s</a></td>' % (self.taxon.id, self.taxon.scientific_name)
        if self.cluster_size == 1:
            ret_string += '</tr></table>'
        else:
            ret_string += '<td><a href="javascript:;" class="cluster-link tip" title="Click to view cluster members" id="%d">(%d)</a></td></tr></table>' % (self.id, self.cluster_size)
        return ret_string

    def get_cluster_pid_link(self):
        return '<a href="javascript:;" id="%d" class="candidate-ortholog-alignment">%.1f</a>' % (self.cluster_representative.id, self.cluster_representative.member_to_query_pid)

    def get_swissprot_link(self):
        if self.cluster_representative.uniprot.in_swissprot_f:
            return '<img src="/static/img/icons/icon_swiss_flag_12.png" />'
        else:
            return ' '

    def get_candidate_ortholog_row_datatable_object(self):
        try:
            return self._get_candidate_ortholog_row_datatable_object
        except:
            pass
        
        gene_name = UniProtGene.objects.filter(uniprot = self.cluster_representative.uniprot)
        try:
            gene_name = gene_name[0].name
        except:
            gene_name = ''

        self._get_candidate_ortholog_row_datatable_object = {
            'DT_RowId': 'co_%d' % self.id,
            '0': self.get_swissprot_link(),
            '1': self.get_uniprot_id_link(),
            '2': self.cluster_representative.uniprot.description,
            '3': gene_name,
            '4': self.get_cluster_taxonomy_link(),
            '5': self.get_cluster_pid_link(),
            '6': '%.1f' % self.cluster_representative.query_coverage,
            '7': '%.1f' % self.cluster_representative.member_coverage
        }
        return self._get_candidate_ortholog_row_datatable_object
    
    class Meta:
        db_table = u'fatcat2_job_member_clusters'

class FatcatJobMembers(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name='members')
    uniprot = models.ForeignKey(UniProt)
    alignment = models.TextField()
    query_coverage = models.FloatField()
    member_coverage = models.FloatField()
    member_to_query_pid = models.FloatField()
    is_query = models.BooleanField()
    is_ortholog = models.BooleanField()
    reason_eliminated = models.TextField()
    cluster = models.ForeignKey(FatcatJobMemberClusters)
    is_cluster_rep = models.BooleanField()
    classification = models.TextField()

    def get_member_header(self):
        if self.is_query:
            return '>lcl|QUERY| ' + self.fatcat_job.fasta_header
        if self.uniprot.in_swissprot_f:
            db = 'sp'
        else:
            db = 'tr'
        header = '>' + db + '|' + self.uniprot.accession + '|' + self.uniprot.uniprot_identifier + ' ' + self.uniprot.description + ' OS=%s' % self.uniprot.taxon.scientific_name
        genes = UniProtGene.objects.filter(uniprot = self)
        if genes:
            header += ' GN=%s' % genes[0].name
        return header

    def get_full_alignment(self):
        sequence = ''
        for (idx, char) in enumerate(self.alignment):
            if ((idx % 60) == 59):
                sequence += '\n'
            sequence += char
        return self.get_member_header() + '\n' + sequence
        
    def get_uniprot_id_link(self):
        return '<a href="http://www.uniprot.org/uniprot/%s" target="_blank">%s</a>' % (self.uniprot.accession, self.uniprot.uniprot_identifier)

    def get_reason_eliminated_link(self):
        if (self.is_ortholog):
            return '<i class="icon-ok tip" title="Meets orthology criteria"></i>'
            #return '<span style="color:green; font-size:15px; font-weight:bold;" class="tip" title="Meets orthology criteria">!</span>'
        else:
            # this is incredibly stupid :'(
            if self.reason_eliminated.startswith('In'):
                return '<i class="icon-ok tip" title="%s"></i>' % self.reason_eliminated
                #return '<span style="color:green; font-size:15px; font-weight:bold;" class="tip" title="%s">!</span>' % self.reason_eliminated
            else:
                return '<i class="icon-remove tip" title="%s"></i>' % self.reason_eliminated
                #return '<span style="color:red; font-size:15px;" class="tip" title="%s">!</span>' % self.reason_eliminated

    def get_member_taxonomy_link(self):
        if self.uniprot.taxon.common_name:
            return '<a href="http://www.uniprot.org/taxonomy/%d" class="tip" title="%s" target="_blank">%s</a>' % (self.uniprot.taxon.id, self.uniprot.taxon.common_name, self.uniprot.taxon.scientific_name)
        else:
            return '<a href="http://www.uniprot.org/taxonomy/%d" target="_blank">%s</a>' % (self.uniprot.taxon.id, self.uniprot.taxon.scientific_name)

    def get_cluster_taxonomy_link(self):
        #if (FatcatJobMember.objects.filter(fatcat_job=self.fatcat_job, cluster_num=self.cluster_num).count() > 1)len(self.fatcat_job.cluster
        if self.uniprot.taxon.common_name:
            return '<a href="http://www.uniprot.org/taxonomy/%d" class="tip" title="%s" target="_blank">%s</a>' % (self.uniprot.taxon.id, self.uniprot.taxon.common_name, self.uniprot.taxon.scientific_name)
        else:
            return '<a href="http://www.uniprot.org/taxonomy/%d" target="_blank">%s</a>' % (self.uniprot.taxon.id, self.uniprot.taxon.scientific_name)
    
    def get_member_pid_link(self):
        return '<a href="javascript:;" id="%d" class="candidate-ortholog-alignment">%.1f</a>' % (self.id, self.member_to_query_pid)

    def get_swissprot_link(self):
        if self.uniprot.in_swissprot_f:
            return '<img src="/static/img/icons/icon_swiss_flag_12.png" />'
        else:
            return ' '

    def get_candidate_ortholog_row_datatable_object(self):
        try:
            return self._get_candidate_ortholog_row_datatable_object
        except:
            pass
        
        gene_name = UniProtGene.objects.filter(uniprot = self.uniprot)
        try:
            gene_name = gene_name[0].name
        except:
            gene_name = ''

        self._get_candidate_ortholog_row_datatable_object = {
            'DT_RowId': 'co_%d' % self.id,
            '0': self.get_swissprot_link(),
            '1': self.get_uniprot_id_link(),
            '2': self.uniprot.description,
            '3': gene_name,
            '4': self.get_member_taxonomy_link(),
            '5': self.get_member_pid_link(),
            '6': '%.1f' % self.query_coverage,
            '7': '%.1f' % self.member_coverage
        }
        return self._get_candidate_ortholog_row_datatable_object
        

    class Meta:
        db_table = u'fatcat2_job_members'

class FatcatJobGOAnnotation(models.Model):
    ''' This is where we store the GO annotations model.  '''
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name='go')
    member = models.ForeignKey(FatcatJobMembers, related_name='go_annotations')
    annotation = models.ForeignKey(Annotation, related_name='member_go_annotation')
   
    def get_annotation_evidence_string(self):
        return '%s: %s' % (self.annotation.evidence.code, self.annotation.evidence.name)

    def get_summary_row_datatable_object(self):
        ''' Table row for summary go annotation for this record. '''
        a = {
            '0': self.member.get_uniprot_id_link(),
            '1': self.member.uniprot.description,
            '2': self.member.get_member_taxonomy_link(),
            '3': self.get_annotation_evidence_string()
        } 
        return a
 
    class Meta:
        db_table = u'fatcat2_job_go_annotation'

''' These are for later.
class FatcatJobTreeNode(models.Model):
    id = models.AutoField(primary_key=True)
    fatcat_job = models.ForeignKey(FatcatJob, related_name='tree_nodes')
    left = models.IntegerField()
    right = models.IntegerField()
    parent_left = models.IntegerField()
    minimum_pairwise_identity_belvu = models.FloatField()
    member = models.ForeignKey(FatcatJobMembers) 
    distance_to_parent = models.FloatField()
    bootstrap_support = models.FloatField()
    sh_test_score = models.FloatField()
    is_query = models.BooleanField()

    class Meta:
        db_table = u'fatcat2_job_tree_node'


class FatcatJobMemberFamilyUnion(models.Model):
    id = models.AutoField(primary_key=True)
    member = models.ForeignKey(FatcatJobMembers)
    member_family = models.ForeignKey(FatcatJobFamily)

    class Meta:
        db_table = u'fatcat2_job_member_family_union'
'''
