from django.db import models
import sys, os, tempfile, pickle, shutil
from pfacts003.intrepid.consts import *
from pfacts003.intrepid.utils import mafft, fast_tree, intrepid, hmmsearch, hmmbuild, \
                                    new_get_fasta_for_accession, get_query_rec, mafft_files, \
                                    get_pdb_file_for_id, get_pdb_fasta_for_id, discern, get_fasta_for_accession, \
                                    get_msa_rec, hmmalign, convert_a2m, add_reference_annotation_line

from pfacts003.phylofacts.models import OhanaJobs, PhyloFactsEmail
from pfacts003.utils.hmm import pfamAscan, hmmscan, parse_domtblout
from bpg.common.BPGPWID import pairwise_identity_belvu as pairwise_identity
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

class CSA(models.Model):
    id = models.IntegerField(primary_key=True)
    pdb_id = models.TextField(db_column='acc')
    chain_id = models.TextField(db_column='chain')
    residue_number = models.IntegerField(db_column='resnum')
    residue = models.TextField()
    function = models.TextField()
    evidence = models.TextField()
    source = models.TextField()
    version = models.TextField()

    class Meta:
        db_table = u'csa'


class QuerySequence(object):
    def __init__(self, p, seq, j, d, dom):
        self.path = p
        self.sequence = seq
        self.job = j
        self.is_domain = d
        self.domain = dom
        self.num_iter = 0
        self.run_discern = False
        self.uniref90_accs = []
        self.prepender = ''
        if self.job.seq_db == UNIPROT_ALL or self.job.seq_db == UNIPROT_ALL_OLD:
            self.prepender = ''
        elif self.job.seq_db == UNIREF100 or self.job.seq_db == UNIREF100_OLD:
            self.prepender = 'UniRef100_'
        elif self.job.seq_db == UNIREF90 or self.job.seq_db == UNIREF90_OLD:
            self.prepender = 'UniRef90_'

    ###################### FILE PATHS ###############################
    def homolog_alignment_path(self, iter=0):
        return os.path.join(self.path, 'homolog_alignment.%d.afa' % (iter))

    def homolog_hmm_path(self, iter=0):
        return os.path.join(self.path, 'homolog_hmm.%d' % (iter))

    def hmmsearch_domain_output_path(self, iter=0):
        return os.path.join(self.path, 'domain.%d.out' % (iter))

    def member_fasta_path(self):
        return os.path.join(self.path, 'member.fasta')

    def percent_id_pickle_path(self, iter=0):
        return os.path.join(self.path, 'percent_ID.%d.pickle' % (iter))
    
    def tree_path(self):
        return os.path.join(self.path, 'member_tree.newick' )

    def tree_alignment_path(self):
        return os.path.join(self.path, 'homolog_alignment_tree.afa')

    def intrepid_alignment_path(self):
        return os.path.join(self.path, 'homolog_alignment_intrepid.afa')

    def intrepid_config_path(self):
        return os.path.join(self.path, 'intrepid.cfg')

    def intrepid_out_path(self):
        return os.path.join(self.path, 'intrepid.out')

    def temp_stockholm_path(self, iter=0):
        return os.path.join(self.path, 'msa.%d.stockholm' % (iter))

    def get_homolog_hmmsearch_options(self):
        return ['-E', str(self.job.homolog_evalue), '--domE', str(self.job.homolog_evalue)]

    def new_get_homologs(self):
        fasta_file = SEQ_DBS[self.job.seq_db]['fasta']
        bdb = SEQ_DBS[self.job.seq_db]['blastable_db']
        if self.is_domain:
            ret_obj = self.domain
            ret_obj.substatus += 'Gathering homologs on %s against %s\n' % (ret_obj.pfam_shortname, fasta_file)
        else:
            ret_obj = self.job
            ret_obj.substatus += 'Gathering homologs on full length protein against %s\n' % (fasta_file)
        # this one is going
        ret_obj.is_running = True
        ret_obj.save()
        # log the start
        num_iter = 0
        divergence = 1
        # write the query to file as homolog_alignment 0
        f = open( self.homolog_alignment_path(num_iter), 'wb')
        f.write(self.sequence)
        f.close()
        # set the current alignment path to be the query file
        curr_ali_path = f.name
        # beginning of the member fasta file
        member_fasta_file = open( self.member_fasta_path(), 'wb')
        member_fasta_file.write(self.sequence+ '\n')
        member_fasta_file.close()

        # flag to terminate if we have gathered no new homologs in this iteration
        terminate = False
        # set of homologs we have already seen
        seen = set()
        return None

    def get_homologs(self):
        ''' this function will get homologs to the query sequence. '''
        # get correct seq_db
        fasta_file = SEQ_DBS[self.job.seq_db]['fasta']
        bdb = SEQ_DBS[self.job.seq_db]['blastable_db']
        if self.is_domain:
            ret_obj = self.domain
            ret_obj.substatus += 'Gathering homologs on %s against %s\n' % (ret_obj.pfam_shortname, fasta_file)
        else:
            ret_obj = self.job
            ret_obj.substatus += 'Gathering homologs on full length protein against %s\n' % (fasta_file)
        # this one is going
        ret_obj.is_running = True
        ret_obj.save()
        # log the start
        num_iter = 0
        divergence = 1
        # write the query to file as homolog_alignment 0
        f = open( self.homolog_alignment_path(num_iter), 'wb')
        f.write(self.sequence)
        f.close()
        # set the current alignment path to be the query file
        curr_ali_path = f.name
        # beginning of the member fasta file
        member_fasta_file = open( self.member_fasta_path(), 'wb')
        member_fasta_file.write(self.sequence+ '\n')
        member_fasta_file.close()

        # flag to terminate if we have gathered no new homologs in this iteration
        terminate = False
        # set of homologs we have already seen
        seen = set()

        # loop if we haven't reached our maximum iterations and if the divergence is still higher than we'd like
        # and if we have anything new to add (otherwise we just do the same hmmsearch over and over again)
        while (num_iter < self.job.maximum_homolog_iterations) and (divergence > self.job.homolog_divergence) and (not terminate): 
            msa = None
            # do the ith iteration
            num_iter += 1
            # log the job
            ret_obj.substatus += 'Starting homolog gathering iteration %d\n' % (num_iter)
            ret_obj.save()
            # build hmm of current iteration's alignment
            hmm_file = open(self.homolog_hmm_path(num_iter),'w')
            add_reference_annotation_line(curr_ali_path, self.temp_stockholm_path(num_iter))
            hmm_file.write(hmmbuild(self.temp_stockholm_path(num_iter), 'ali%d' % (num_iter), ['--hand']))
            hmm_file.close()
            # do the hmmsearch            
            hmmsearch(hmm_file.name, fasta_file, self.hmmsearch_domain_output_path(num_iter), 
                            self.get_homolog_hmmsearch_options())
            # parse the output file and sort according to bit score
            sorted_rows = sorted(parse_domtblout(open(self.hmmsearch_domain_output_path(num_iter),'rb').read()), 
                                    key=lambda r: r.full_score, reverse=True)

            final_rows = []

            for row in sorted_rows:
                # need different splitting protocol because fasta header line is different
                # for different versions of uniprot.  this could be done way better
                if self.job.seq_db == UNIPROT_ALL or self.job.seq_db == UNIPROT_ALL_OLD:
                    try:
                        acc = row.target_name.split('|')[1]
                    except:
                        continue
                else:
                    try:
                        acc = row.target_name
                    except:
                        continue
                if (acc not in seen) and float(1.0*(int(row.hmm_to)-int(row.hmm_from)))/int(row.query_length) >= self.job.homolog_coverage_condition:
                    seen.add(acc)
                    final_rows.append(row)
                    if 'UniRef' in acc:
                        new_acc = acc.split('_')[1]
                    else:
                        new_acc = acc
                    mf = get_fasta_for_accession('UniRef90_' + new_acc, SEQ_DBS[UNIREF90_OLD]['blastable_db'])
                    if mf.startswith('>'):
                        self.uniref90_accs.append(new_acc)
            if not final_rows:
                # break out of here, we have nothing more to do, set the flag to terminate.
                ret_obj.substatus += 'No new homologs were found.  Exiting.\n'
                ret_obj.save()
                terminate = True
            # make the fasta file of the query plus the homologs    
            qseq = ''.join(self.sequence.split('\n')[1:])
            # append to that vial
            member_fasta_file = open(self.member_fasta_path(), 'a')
            for row in final_rows:
                if self.job.seq_db == UNIPROT_ALL or self.job.seq_db == UNIPROT_ALL_OLD:
                    acc = row.target_name.split('|')[1]
                else:
                    acc = row.target_name
                if self.job.use_old_blasttools:
                    mf = get_fasta_for_accession(acc, bdb)
                else:
                    mf = new_get_fasta_for_accession(acc, bdb)
                hdr = mf.split('\n')[0]
                seq = ''.join(mf.split('\n')[1:])
                # if its a domain, get the aligned section
                if self.is_domain:
                    new_header = hdr + ' %d,%d' % (row.ali_from, row.ali_to)
                    new_sequence = seq[row.ali_from-1:row.ali_to]
                else:
                    new_header = hdr
                    new_sequence = seq
                if (qseq != new_sequence):
                    member_fasta_file.write(new_header+'\n'+new_sequence+'\n')
            member_fasta_file.close()
            final_rows = None
            
            # align the homologs with hmmalign and write it to file
            f = open(self.homolog_alignment_path(num_iter), 'wb')           
            f.write(convert_a2m(hmmalign(member_fasta_file.name, hmm_file.name, options=['--informat', 'FASTA', '--outformat', 'afa'])))
            f.close()
            # set the current alignment path
            curr_ali_path = self.homolog_alignment_path(num_iter)
            # load the msa into a biopython alignment
            msa = AlignIO.read(open(self.homolog_alignment_path(num_iter),'rb'), 'fasta')
            # query record
            qrec = get_query_rec(msa)
            # save the percent ID structure
            percent_ID = {}
            # now find the minimum pairwise identity
            for rec in msa:
                if rec != qrec:
                    pwid = pairwise_identity(str(rec.seq), str(qrec.seq))
                    # to pickle the data structure
                    percent_ID[rec.description] = pwid
                    if pwid < divergence:
                        divergence = pwid
            # pickle the percent ID so we don't have to recalculate this.
            pickle.dump(percent_ID, open(self.percent_id_pickle_path(num_iter), 'wb'))
            ret_obj.substatus += 'At the end of iteration %d: %d homologs were found with a maximum divergence of %.2f\n' % (num_iter, (len(msa)-1), divergence)
            ret_obj.save()
        self.num_iter = num_iter
        ret_obj.substatus += 'Ending homolog gathering in %d iterations with %d homologs and a maximum divergence from the query of %.2f\n' % (num_iter, (len(msa)-1), divergence)
        ret_obj.is_running = False
        ret_obj.save()
        return (num_iter, divergence, (len(msa)-1), ret_obj)

    def process_msa(self):
        # this will process the final homolog gathering msa.
        f = open(self.homolog_alignment_path(self.num_iter), 'rb')
        msa = AlignIO.read(f, 'fasta')
        '''if self.job.is_discern_job:
            qrec = get_query_rec(msa, qid=self.job.pdb_id())
        else:
            qrec = get_query_rec(msa)'''
        qrec = get_query_rec(msa)
        intrepid_msa = []
        tree_msa = []

        # if its a domain, set ret_obj to the domain, otherwise to the job....
        if self.is_domain:
            ret_obj = self.domain
        else:
            ret_obj = self.job
    
        # first remove any exact duplicates of the seed, screw up tree builders
        #msa = MultipleSeqAlignment([rec for rec in msa if ((rec == qrec) or (str(rec.seq) != str(qrec.seq)))]) 

        num_rows = len(msa)

        # next loop through the columns and do the filtering magic
        for column in range(len(str(qrec.seq))):
            gappiness = 1.0*msa[:,column].count('-')/num_rows
            # first see if the column is sufficiently ungappy
            if (gappiness < self.job.gappy_column_condition):
                if tree_msa:
                    tree_msa = tree_msa + msa[:,column:(column+1)]
                else:
                    tree_msa = msa[:,column:(column+1)]
            # next see if the query is gapped here, this msa goes to intrepid, as the variable name suggests
            if (str(qrec.seq)[column] != '-'):
                if intrepid_msa:
                    intrepid_msa = intrepid_msa + msa[:,column:(column+1)]
                else:
                    intrepid_msa = msa[:,column:(column+1)]
        qrec = get_query_rec(intrepid_msa)
        trec = get_query_rec(tree_msa)
        new_intrepid_msa = [qrec]
        new_tree_msa = [trec]
        removed = 0
        percent_ID = {}
        removed_list = ''
        lowest_ID = 1000
        
        # remove the entries that don't pass the coverage condition forcibly.
        for rec in intrepid_msa:
            if rec != qrec:
                if (float(1.0*sum([1 for x in str(rec.seq) if x.isupper()]))/(len(intrepid_msa[0])) >= \
                    self.job.homolog_coverage_condition):
                    new_intrepid_msa.append(rec)
                    # get the corresponding record in the tree msa
                    for r in tree_msa:
                        if r.description == rec.description:
                            new_tree_msa.append(r)
                            break
                    pwid = pairwise_identity(str(r.seq), str(trec.seq))
                    # to pickle the data structure
                    percent_ID[r.description] = pwid
                    if pwid < lowest_ID:
                        lowest_ID = pwid
                else:
                    removed += 1
                    removed_list += '%s, ' % (rec.id)
        # pickle the percent ID so we don't have to recalculate this.
        intrepid_msa = MultipleSeqAlignment(new_intrepid_msa)
        tree_msa = MultipleSeqAlignment(new_tree_msa)
        #ret_obj.substatus += 'Removed %d (%.2f %%) sequences not meeting coverage criteria.  %d sequences remain.\n' % (removed, (100.0*removed)/(removed + len(intrepid_msa)), len(intrepid_msa)-1)
        ret_obj.num_homologs -= removed
        ret_obj.maximum_divergence = lowest_ID
        #ret_obj.substatus += 'Removed accessions: %s\n' % removed_list
        ret_obj.substatus += '%s sequences remain with a percent ID between most remote homolog and query of %.2f %%\n' % (ret_obj.num_homologs, 100.0*ret_obj.maximum_divergence)
        ret_obj.save()
        if (len(intrepid_msa) - 1) > self.job.max_homologs:
            ret_obj.substatus += 'Too many homologs.  Starting to cull remaining.\n'
            ret_obj.save()
            qrec = get_query_rec(intrepid_msa)
            trec = get_query_rec(tree_msa)
            new_intrepid_msa = [qrec]
            new_tree_msa = [trec]
            count = 0
            for acc in self.uniref90_accs:
                nim = get_msa_rec(intrepid_msa, self.prepender + acc)
                ntm = get_msa_rec(tree_msa, self.prepender + acc) 
                if nim and ntm:
                    new_intrepid_msa.append(nim)
                    new_tree_msa.append(ntm)
                    count += 1
                if count == self.job.max_homologs:
                    break
            intrepid_msa = MultipleSeqAlignment(new_intrepid_msa)
            tree_msa = MultipleSeqAlignment(new_tree_msa)
            trec = get_query_rec(tree_msa)
            percent_ID = {}
            divergence = 1
            # get new min pid
            for rec in tree_msa:
                if rec != trec:
                    pwid = pairwise_identity(str(rec.seq), str(trec.seq))
                    percent_ID[rec.description] = pwid
                    if pwid < divergence:
                        divergence = pwid
            ret_obj.num_homologs = len(tree_msa) - 1
            ret_obj.maximum_divergence = divergence
            ret_obj.substatus += 'Sequence cull done.  %d homologs remain with a percent ID between query and most remote homolog of %.2f\n' % (len(tree_msa), 100.0*divergence)
            ret_obj.save()
        pickle.dump(percent_ID, open(self.percent_id_pickle_path(ret_obj.num_iterations), 'wb'))
        g = open(self.intrepid_alignment_path(),'wb')
        h = open(self.tree_alignment_path(), 'wb')
        AlignIO.write(intrepid_msa, g, 'fasta')
        AlignIO.write(tree_msa, h, 'fasta')
        # log the information about the number of sequences removed
        f.close()
        g.close()
        h.close()
        return (g.name, h.name)
    
    def build_tree(self):
        f = open(self.tree_path(), 'wb')
        f.write(fast_tree(self.tree_alignment_path()))
        f.close()
        return None

    def run_intrepid(self):
        # make the config file for intrepid
        icf = open(self.intrepid_config_path(), 'wb')
        icf.write("msa_file %s\ntree_file %s\nsequence_id QUERY\n" % (self.intrepid_alignment_path(), self.tree_path()))
        icf.close()
        # change the dir for intrepid
        os.chdir(self.path)
        # then run it
        intrepid(self.intrepid_config_path(), self.intrepid_out_path())
        return None

class IntrepidStatus(models.Model):
    id = models.AutoField(primary_key=True)
    pipeline_stage = models.IntegerField()
    total_pipeline_stages = models.IntegerField()
    status = models.TextField()

    def status_percent(self):
        return float(1.0*self.pipeline_stage)/self.total_pipeline_stages

    class Meta:
        db_table = u'intrepid_status'

class IntrepidJob(models.Model):
    id = models.AutoField(primary_key=True)
    # fasta stuff
    fasta_header = models.TextField()
    fasta_sequence = models.TextField()
    # pdb id - holds the pdb id, or the path to the uploaded pdb file.  If this is none
    # it is just an intrepid job, if this exists it is a discern job also
    pdb_id = models.TextField()
    # this is the name of this job, if its an intrepid job, this is the whole protein, if its discern
    # this is the pdb_id + the chain information
    internal_name = models.TextField()

    # job logging stuff
    ohana_job = models.ForeignKey(OhanaJobs, related_name = 'fsp_job')
    # status counter
    status = models.ForeignKey(IntrepidStatus, related_name = 'fsp_status')

    ###########################
    # parameters for this job
    ###########################
    # flag if we should save this job
    save_job = models.BooleanField()
    # flag for if the job is done
    is_done = models.BooleanField()
    # pbs job id
    pbs_job_id = models.TextField()
    # path where the job files can be found
    job_path = models.TextField()

    # homolog E-value
    homolog_evalue = models.FloatField()

    # sequence database selector
    seq_db = models.IntegerField()

    # maximum number of iterations for homolog gathering
    maximum_homolog_iterations = models.IntegerField()
    # divergence necessary to stop homolog gathering
    homolog_divergence = models.FloatField()
    # remove columns with greater than this percentage gaps
    gappy_column_condition = models.FloatField()
    # minimum pfam domain length for analysis
    minimum_pfam_length = models.IntegerField()    
    # is this a discern job?
    is_discern_job = models.BooleanField()
    # the minimum coverage required for inclusion into homologs
    homolog_coverage_condition = models.FloatField()
    # should we use the blastdb created by the older blast tools?  formatdb and fastacmd
    use_old_blasttools = models.BooleanField()
    # maximum number of homologs allowed
    max_homologs = models.IntegerField()
    # should we save the program outputs to a file?
    save_program_outputs = models.BooleanField()

    ###########################
    # Job Stats
    ###########################
    # num of homologs found
    num_homologs = models.IntegerField()
    # number of iterations required
    num_iterations = models.IntegerField()
    # maximum divergence
    maximum_divergence = models.FloatField()
    # substatus for full length
    substatus = models.TextField()
    # is this job currently running ?
    is_running = models.BooleanField()
    # what is the user supplied name for this job?
    job_name = models.TextField()
    # holds the alignment of the sequence obtained from pdb fasta file
    pdb_fasta_alignment_line = models.TextField()
    # holds the alignment of the sequence obtained from the pdb structure file
    pdb_structure_alignment_line = models.TextField()
    # holds the offset between the pdb structure sequence numbering and the fasta numbering
    pdb_structure_fasta_offset = models.IntegerField()
    # does this structure contain internal gaps?
    structure_contains_internal_gaps = models.BooleanField()
    # does this job have an uploaded model?
    has_uploaded_model = models.BooleanField()
    # what is the comparative model string?
    comparative_model_string = models.TextField()

    ###########################
    # email stuff
    ###########################
    user_email = models.TextField()

    ###########################
    # is this a development job?
    ###########################
    development_job = models.BooleanField(default=False)
    
    def __init__(self, *args, **kwargs):
        models.Model.__init__(self, *args, **kwargs)
        # class variable for the queries that we run
        self.queries = []
        # try to get the pdb file path, set the intrepid job variable
        if self.pdb_id:
            '''try:
                f = open(self.pdb_id,'rb')
                self.pdb_file_path = f.name
                f.close()
            except:
                self.pdb_file_path = open(get_pdb_file_for_id(id))
            '''
            self.is_intrepid_job = False            
        else:
            self.is_intrepid_job = True

    @property
    def fasta(self):
        if self.fasta_header.startswith('>'):
            return self.fasta_header + '\n' + self.fasta_sequence
        else:
            return '>' + self.fasta_header + '\n' + self.fasta_sequence

    @property
    def name(self):
        # this is the name of the left tab for the job
        if self.is_discern_job:
            return self.pdb_id[:4]
        else:
            return 'Whole Protein'

    @property
    def chain(self):
        # this is the name of the chain if its a discern job
        if self.is_discern_job:
            return self.pdb_id[4:]
        else:
            return None

    ##############################
    # FILE PATHS
    ##############################
    @property
    def homolog_alignment_tree_path(self):
        return os.path.join(self.job_path, 'homolog_alignment_tree.afa')

    @property
    def homolog_alignment_intrepid_path(self):
        return os.path.join(self.job_path, 'homolog_alignment_intrepid.afa')
    @property
    def tree_path(self):
        return os.path.join(self.job_path, 'member_tree.newick')

    @property
    def intrepid_config_path(self):
        return os.path.join(self.job_path, 'intrepid.cfg')

    @property
    def intrepid_out_path(self):
        return os.path.join(self.job_path, 'intrepid.out')

    @property
    def intrepid_rank_file_path(self):
        return os.path.join(self.job_path, 'output.rank')

    @property
    def intrepid_aux_file_path(self):
        return os.path.join(self.job_path, 'output.aux')
   
    @property
    def discern_output_file_path(self):
        return os.path.join(self.job_path, self.name, 'INTREPID', '%s.discern' % self.name)
    
    @property
    def discern_centrality_file_path(self):
        return os.path.join(self.job_path, self.name, 'INTREPID', '%s.central' % self.name)

    @property
    def discern_zipfile_path(self):
        return os.path.join(self.job_path, 'Job%d_discern.gzip' % (self.id)) 
    
    @property
    def comparative_model_path(self):
        return os.path.join(self.job_path, '%s' % self.name, 'INTREPID','%s.pdb' % self.name)

    @property
    def original_pdb_file_path(self):
        return os.path.join(self.job_path, '%s' % self.name, 'INTREPID', '%s.original.pdb' % self.name)

    @property
    def program_outputs_path(self):
        return os.path.join(self.job_path, 'program_outputs.txt')

    @property
    def most_remote_homolog_info(self):
        f = open(self.percent_ID_pickle_path, 'rb')
        pids = pickle.load(f)
        f.close()
        return sorted(pids.items(), key=lambda x: x[1])[0]        

    @property
    def percent_ID_pickle_path(self):
        return os.path.join(self.job_path, 'percent_ID.%d.pickle' % (self.num_iterations))

    def pretty_sequence(self, length=60):
        # pretty prints the sequence
        index = 0
        retstr = ''
        for char in self.fasta_sequence:
            if index == length:
                retstr += '\n'
                index = 0
            retstr += char
            index += 1
        return retstr

    def status_object(self):
        ret_obj = {}
        # returns the status of the job for the progress page, lots of redundant data sent here
        ret_obj['full_length'] = self.substatus
        for domain in self.pfams_analyzed():
            ret_obj['D_%d' % (domain.id)] = {
                'title': domain.pfam_shortname + ' domain',
                'status': domain.substatus
            }
        ret_obj['id'] = self.id
        ret_obj['status'] = self.status.status
        ret_obj['stage'] = self.status.pipeline_stage
        ret_obj['total_stages'] = self.status.total_pipeline_stages
        return ret_obj

    def get_csa_residues(self):
        # return a comma separated string with the csa residues 
        return ','.join([str(r.residue_number) for r in CSA.objects.filter(pdb_id = self.pdb_id[0:4], chain_id = self.pdb_id[4:])])
            
    def get_intrepid_ranks(self):
        # return a vector of structure positions in order of intrepid rank
        ret= []
        for res in self.residues.all().order_by('intrepid_rank'):
            ret.append([res.structure_position, res.intrepid_cons_js_score])
        return ret

    def get_intrepid_scores(self):
        # return a vector of structure positions in order of intrepid rank
        ret= []
        for res in self.residues.all().order_by('-intrepid_score'):
            ret.append(res.intrepid_cons_js_score)
        return ret

    def get_discern_ranks(self):
        # return a vector of structure positions in order of discern rank
        ret = []
        for res in self.residues.all().order_by('discern_rank'):
            ret.append([res.structure_position, res.discern_score])
        return ret
    
    def get_discern_scores(self):
        # return a vector of structure positions in order of discern rank
        ret = []
        for res in self.residues.all().order_by('-discern_score'):
            ret.append(res.discern_score)
        return ret

    def intrepid_scores(self):
        # memoize reading these scores...
        try:
            return self._intrepid_scores
        except:
            pass
        f = open(self.intrepid_aux_file_path, 'r')
        scores = []
        for line in f.readlines()[1:]:
            fields = line.split('|')
            scores.append(( "%s (%s)" % (str(fields[3]), AA_DICT[str(fields[3])]['3']), 
                    AA_DICT[str(fields[3])]['full'], int(fields[0]), float(fields[6]), float(fields[9]) ))
        f.close()
        self._intrepid_scores = sorted(scores, key=lambda x: x[3],reverse=True)
        return self._intrepid_scores

    def intrepid_chart_series(self):
        return [ [x[2], x[3]] for x in self.intrepid_scores() ]

    def discern_scores(self):
        try:
            self._discern_scores
        except:
            pass
        f = open(self.discern_output_file_path, 'rb')
        self._discern_scores = sorted([ ( int(line.split()[0]), float(line.split()[1]), \
            self.fasta_sequence[int(line.split()[0]) - 1], AA_DICT[self.fasta_sequence[int(line.split()[0]) - 1]]['3'], \
            AA_DICT[self.fasta_sequence[int(line.split()[0])-1]]['full']) for line in f.readlines() \
            if int(line.split()[0]) <= len(self.fasta_sequence)], key = lambda x: x[1], reverse=True)
        f.close()
        return self._discern_scores

    def scan_pfams(self):
        self.job_path = tempfile.mkdtemp(dir='/clusterfs/ohana/software/webserver/incoming/')
        self.save()
        # make the necessary directories
        seed_dir = os.path.join(self.job_path, self.name)
        os.mkdir(seed_dir)
        intrepid_dir = os.path.join(seed_dir, 'INTREPID')
        os.mkdir(intrepid_dir)
        # write the comparative model to string
        if self.has_uploaded_model:
            f = open(self.comparative_model_path, 'wb')
            f.write(self.comparative_model_string)
            f.close()
            self.comparative_model_string = ''
            self.save()
        pfam_hits = pfamAscan(self.fasta, options=['-E','1e-3'])
        for pfam in pfam_hits:
            intrepid_job_pfam = self.pfams.create(
                    pfam_description = pfam.description,
                    pfam_accession = pfam.target_accession,
                    pfam_shortname = pfam.target_name,
                    ali_from = pfam.ali_from,
                    ali_to = pfam.ali_to,
                    hmm_from = pfam.hmm_from,
                    hmm_to = pfam.hmm_to,
                    c_evalue = pfam.this_domain_c_value,
                    i_evalue = pfam.this_domain_i_value,
                    analyzed = False,
                    num_homologs = 0,
                    num_iterations = 0,
                    maximum_divergence = 0.0,
                    substatus = 'Initialized.\n',
                    is_running = False
                )
            if ((pfam.ali_to-pfam.ali_from) >= self.minimum_pfam_length) and not self.is_discern_job:
                # Make the temporary dir for the domains
                os.mkdir(os.path.join(self.job_path, str(intrepid_job_pfam.id)))
                intrepid_job_pfam.analyzed = True
                intrepid_job_pfam.substatus += 'Minimum length of Pfam domain for separate analysis: %d\n' % (self.minimum_pfam_length)
                intrepid_job_pfam.substatus += 'ID necessary between query and most remote homolog: %d%%\n' % (int(100 * self.homolog_divergence))
                intrepid_job_pfam.substatus += 'Minimum query (hmm) coverage for inclusion into homolog set: %d%%\n' % (int(100 * self.homolog_coverage_condition))
                intrepid_job_pfam.substatus += 'Maximum iterations of homolog gathering: %d\n' % (self.maximum_homolog_iterations)
                intrepid_job_pfam.substatus += 'Minimum column coverage for inclusion into tree MSA: %d%%\n'% (int(100 * self.gappy_column_condition))
                intrepid_job_pfam.substatus += 'Maximum number of homologs: %d\n' % (self.max_homologs)
                intrepid_job_pfam.substatus += 'Sequence database to gather homologs from: %s\n' % (SEQ_DBS[self.seq_db]['name'])
                    
                intrepid_job_pfam.save()
        return None

    def make_query_sequences(self):
        # This will populate the self.queries variable with an array corresponding to each query that needs to be run.
        # first make query sequence on the full length protein
        self.queries = [QuerySequence(p = self.job_path, 
                            seq='>lcl|QUERY| %s\n%s' % (self.fasta_header.strip('>'), self.fasta_sequence),
                            j = self, d=False,dom=None)]

        # now make query sequences for the domains
        for domain in self.pfams_analyzed():
            self.queries.append( QuerySequence(p = domain.domain_path, 
                                    seq = '>lcl|QUERY| %s %d,%d\n%s' % (self.fasta_header.strip('>'), 
                                    domain.ali_from, domain.ali_to, 
                                    self.fasta_sequence[domain.ali_from-1:domain.ali_to]),
                                    j = self, d=True, dom=domain) )
        return self.queries
    
    def gather_homologs(self):
        # make the query sequences
        self.make_query_sequences()
        # run get homologs on each query sequence
        for s in self.queries:
            (ni,md,nh,obj) = s.get_homologs()
            obj.num_homologs = nh
            obj.num_iterations = ni
            obj.maximum_divergence = md
            obj.save()
        return None

    def pfams_analyzed(self):
        return self.pfams.filter(analyzed=True)
    
    def ordered_pfams(self):
        # returns the pfam domains ordered by the alignment start
        return self.pfams.filter(analyzed=True).order_by('ali_from')
    
    def mask_msa(self):
        # creates the two alignments for each of the QuerySequence objects
        for s in self.queries:
            s.process_msa()
        return None
        
    def build_tree(self):
        for s in self.queries:
            s.build_tree()
        return None        

    def run_intrepid(self):
        for s in self.queries:
            s.run_intrepid()
        if self.is_discern_job:
            self.populate_intrepid_scores()
        return None

    def populate_intrepid_scores(self):
        f = open(self.intrepid_aux_file_path, 'r')
        scores = []
        for (pos, line) in enumerate(f.readlines()[1:]):
            fields = line.split('|')
            this_residue = self.residues.get(fasta_position=int(fields[0]))
            this_residue.intrepid_cons_js_score = float(fields[6])
            this_residue.intrepid_global_js_score = float(fields[9])
            this_residue.save()
        f.close()
        # now put in the rank
        for (rank, res) in enumerate(self.residues.all().order_by('-intrepid_cons_js_score')):
            res.intrepid_rank = rank + 1
            res.save()
        return None            

    def run_discern(self):
        # create discern directory structure
        intrepid_dir = os.path.join(self.job_path, self.name, 'INTREPID')
        # copy the intrepid score into that directory
        shutil.copy(self.intrepid_aux_file_path, intrepid_dir)
        # run discern script
        discern(self.name, self.pdb_id, self.job_path, self.pdb_structure_fasta_offset, self.has_uploaded_model)
        return None

    def populate_discern_scores(self):
        # first populate the discern score and centrality
        g = open(self.discern_output_file_path, 'rb')
        f = open(self.discern_centrality_file_path,'rb')
        # get the scores
        h = g.readlines()
        i = f.readlines()
        g.close()
        f.close()
        for (num, line) in enumerate(h):
            spos = int(line.split()[0])
            out_score = float(line.split()[1])
            out_centrality = float(i[num].split()[1])
            try:
                this_residue = self.residues.get(structure_position=spos+self.pdb_structure_fasta_offset)
                this_residue.discern_centrality_score = out_centrality
                this_residue.discern_score = out_score
                this_residue.save()
            except:
                pass
        # now populate the discern rank
        for (rank, residue) in enumerate(self.residues.all().order_by('-discern_score')):
            residue.discern_rank = rank + 1
            residue.save()
        return None
        
    def make_directories_readable(self):
        os.chmod( self.job_path, 0777 ) 
        for r,d,f in os.walk(self.job_path):
            for dirs in d:
                os.chmod( os.path.join(r,dirs), 0777 )
            for files in f:
                os.chmod( os.path.join(r,files), 0777 )
        return None

    def send_email(self):
        # sends email
        if self.user_email:
            sub = 'Functional Site Prediction Job %d complete' % (self.id)
            text_body = '''Functional Site Prediction Job %d has completed.
                \n\nhttp://phylogenomics.berkeley.edu/intrepid/%d/
                \n
                \n\n%s
                \n\n''' % (self.id, self.id, self.fasta_header)
            html_body = '''<h2>Functional Site Prediction Job %d has completed</h2>.
                \n\n<p><a href="http://phylogenomics.berkeley.edu/intrepid/%d/">click here to view the results</a></p>
                \n\n
                \n\n<p>%s</p>\n
                \n\n''' % (self.id, self.id, self.fasta_header)
            PhyloFactsEmail.objects.create(
                recipient = self.user_email,
                text_body = text_body,
                html_body = html_body,
                subject = sub
            )
        return None        

    def get_substatus_text(self):
        f = open(self.program_outputs_path, 'rb')
        retstr = f.read()
        f.close()
        return retstr

    def log_job_completion(self):
        # convert tree to phyloxml ?
        # write the substatuses to file
        f = open(self.program_outputs_path,'wb')
        f.write(self.substatus+'\n')
        self.substatus = ''
        self.save()
        for d in self.pfams_analyzed():
            f.write(d.substatus+'\n')
            d.substatus = ''
            d.save()
        f.close()
        self.make_directories_readable()
        self.send_email()
        self.is_done = True
        self.status_id = JOB_DONE
        self.save()
        return None

    class Meta:
        db_table = u'intrepid_job'

class IntrepidJobPfamDomain(models.Model):
    id = models.AutoField(primary_key=True)
    intrepid_job = models.ForeignKey(IntrepidJob, related_name="pfams")
    pfam_description = models.TextField()
    pfam_accession = models.TextField()
    pfam_shortname = models.TextField()
    ali_from = models.IntegerField()
    ali_to = models.IntegerField()
    hmm_from = models.IntegerField()
    hmm_to = models.IntegerField()
    c_evalue = models.FloatField()
    i_evalue = models.FloatField()
    
    # should this domain be analyzed (was it analyzed) ?
    analyzed = models.BooleanField()

    ###########################
    # Job Stats
    ###########################
    # num of homologs found
    num_homologs = models.IntegerField()
    # number of iterations required
    num_iterations = models.IntegerField()
    # maximum divergence
    maximum_divergence = models.FloatField()
    # substatus string
    substatus = models.TextField()
    # is this domain currently running?
    is_running = models.BooleanField()

    @property
    def domain_path(self):
    # return the path where the domain files are located
        return os.path.join(self.intrepid_job.job_path, str(self.id))

    ###########################################
    # FILE PATHS
    ###########################################    
    @property
    def homolog_alignment_tree_path(self):
        return os.path.join(self.domain_path, 'homolog_alignment_tree.afa')

    @property
    def homolog_alignment_intrepid_path(self):
        return os.path.join(self.domain_path, 'homolog_alignment_intrepid.afa')

    @property
    def tree_path(self):
        return os.path.join(self.domain_path, 'member_tree.newick')

    @property
    def intrepid_config_path(self):
        return os.path.join(self.domain_path, 'intrepid.cfg')
    
    @property
    def intrepid_out_path(self):
        return os.path.join(self.domain_path, 'intrepid.out')

    @property
    def intrepid_rank_file_path(self):
        return os.path.join(self.domain_path, 'output.rank')

    @property
    def intrepid_aux_file_path(self):
        return os.path.join(self.domain_path, 'output.aux')
    
    @property
    def most_remote_homolog_info(self):
        f = open(self.percent_ID_pickle_path, 'rb')
        pids = pickle.load(f)
        f.close()
        return sorted(pids.items(), key=lambda x: x[1])[0]        

    @property
    def percent_ID_pickle_path(self):
        return os.path.join(self.domain_path, 'percent_ID.%d.pickle' % (self.num_iterations))
    
    @property
    def html_friendly_accession(self):
        return self.pfam_accession.split('.')[0]
    
    def intrepid_scores(self):
        try:
            return self._intrepid_scores
        except:
            pass
        f = open(self.intrepid_aux_file_path, 'r')
        scores = []
        for line in f.readlines()[1:]:
            fields = line.split('|')
            scores.append(( "%s (%s)" % (str(fields[3]), AA_DICT[str(fields[3])]['3']), 
                    AA_DICT[str(fields[3])]['full'], int(fields[0]), float(fields[6]), 
                    float(fields[9]), int(fields[0])+self.ali_from-1 ))
        f.close()
        self._intrepid_scores = sorted(scores, key=lambda x: x[3],reverse=True)
        return self._intrepid_scores

    def intrepid_chart_series(self):
        return [ [x[2], x[3]] for x in self.intrepid_scores() ]

    class Meta:
        db_table = u'intrepid_job_pfam_domain'

class IntrepidJobResidues(models.Model):
    id = models.AutoField(primary_key=True)
    intrepid_job = models.ForeignKey(IntrepidJob, related_name="residues")
    fasta_residue = models.TextField()
    fasta_residue_three_letter = models.TextField()
    structure_residue = models.TextField()
    structure_residue_three_letter = models.TextField()
    fasta_position = models.IntegerField()
    structure_position = models.IntegerField()
    # 
    discern_rank = models.IntegerField()
    #
    discern_score = models.FloatField()
    #
    discern_centrality_score = models.FloatField()
    #
    intrepid_rank = models.IntegerField()
    #
    intrepid_cons_js_score = models.FloatField()
    #
    intrepid_global_js_score = models.FloatField()
    global_conservation = models.FloatField()
    relative_solvent_accessibility = models.FloatField()
    absolute_solvent_accessibility = models.FloatField()

    def support_html(self):
        if not self.intrepid_job.pdb_id:
            return '&nbsp;'
        csa_res = CSA.objects.filter(pdb_id = self.intrepid_job.pdb_id[0:4], chain_id = self.intrepid_job.pdb_id[4:], residue_number = self.structure_position)
        if csa_res:
            return '<a style="color:red;font-size:10px" class="tip" title="Direct Catalytic Site Atlas support" target="_blank" href="http://www.ebi.ac.uk/thornton-srv/databases/CSA/SearchResults.php?PDBID=%s&SUBMIT_PDB=SEARCH+CSA">CSA</a>' % (self.intrepid_job.pdb_id[0:4])
        else:
            return '&nbsp;' 

    class Meta:
        db_table = u'intrepid_job_residues'

class PDB(models.Model):
    id = models.AutoField(primary_key = True)
    pdb_id = models.TextField()
    chain_id = models.TextField()
    type = models.TextField()
    contains_numbering_discontinuities = models.BooleanField()
    sequence_structure_contains_gaps = models.BooleanField()
    sequence_fasta_contains_gaps = models.BooleanField()
    structure_length = models.IntegerField()
    fasta_length = models.IntegerField()
    structure_fasta_offset = models.IntegerField()
    structure_sequence_alignment_line = models.TextField()
    fasta_sequence_alignment_line = models.TextField()
    pdb_download_date = models.TextField()
    num_models = models.IntegerField()
    model_number = models.IntegerField()
    biopython_parse_error = models.BooleanField()

    class Meta:
        db_table = u'new_pdb'

'''
class IntrepidJobStructures(models.Model):
    id = models.AutoField(primary_key=True)
    intrepid_job = models.ForeignKey(IntrepidJob, related_name="structures")
    hmm_from = models.IntegerField()
    hmm_to = models.IntegerField()
    ali_from = models.IntegerField()
    ali_to = models.IntegerField()
    c_evalue = models.FloatField()
    i_evalue = models.FloatField()

    class Meta:
        db_table = 'intrepid_job_structures'
'''
