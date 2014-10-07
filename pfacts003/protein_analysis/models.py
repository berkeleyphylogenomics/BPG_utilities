from django.db import models
from pfacts003.phylofacts.models import OhanaJobs
from pfacts003.utils.hmm import pfamAscan

class ProteinDomainStatus(models.Model):
    id = models.IntegerField(primary_key=True)
    status = models.TextField()

    class Meta:
        db_table = u'protein_domain_status'

class ProteinAnalysisJob(models.Model):
    id = models.AutoField(primary_key=True)
    # fasta stuff
    fasta_header = models.TextField()
    fasta_sequence = models.TextField()
    # job logging stuff
    ohana_job = models.ForeignKey(OhanaJobs, related_name = 'pa_job')
    # status counter
    status = models.IntegerField()    

    ###########################
    # parameters for this job
    ###########################
    
    # minimum aligned length for us to consider a domain
    minimum_domain_length = models.IntegerField()
    # domains we should try to find in the query
    # this is a space separated string containing the names of the domain methods
    # we should try to find in the query.  this will be parsed later.
    domain_methods = models.TextField()
    
    ###########################
    # end parameters
    ###########################
    
    # flag if we should save this job
    save_job = models.BooleanField()
    # flag if the domains of the query have been scanned already
    domains_scanned = models.BooleanField()
    # pbs job id
    pbs_job_id = models.TextField()
    # path where the job files can be found
    job_path = models.TextField()

    @property
    def fasta(self):
        if self.fasta_header.startswith('>'):
            return self.fasta_header + '\n' + self.fasta_sequence
        else:
            return '>' + self.fasta_header + '\n' + self.fasta_sequence

    def populate_pfam_domains(self, method='pfam', full_method_name='Pfam domain'):
        # this function actually scans the query versus pfam hmms and stores the significant
        # non overlapping domains
        pfam_hits = pfamAscan(self.fasta, options=['-E', '1e-3'])
        for pfam in pfam_hits:
            if (pfam.ali_to - pfam.ali_from) >= self.minimum_domain_length:
                self.domains.create(
                    type = method,
                    full_type_name = full_method_name,
                    description = pfam.description,
                    accession = pfam.target_accession,
                    shortname = pfam.target_name,
                    ali_from = pfam.ali_from,
                    ali_to = pfam.ali_to,
                    hmm_from = pfam.hmm_from,
                    hmm_to = pfam.hmm_to,
                    c_evalue = pfam.this_domain_c_value,
                    i_evalue = pfam.this_domain_i_value,
                    homologs_gathered = False,
                    intrepid_done = False,
                    discern_done = False,
                    status_id = 0,
                    currently_running_job = 'None'   
             )
        return None

    def find_query_domains(self):
        ''' this function will try to find domains in the query to build trees on. '''
        ###################################################################################
        # THIS FUNCTION SHOULD ONLY BE CALLED ONCE AS IT POPULATES THE DOMAINS RECORDS
        ###################################################################################
        if self.domains_scanned:
            return None

        # first, go through all of the domain methods we are asked to and make the domain records
        for method in self.domain_methods.strip().split():
            if method == 'pfam':
                self.populate_pfam_domains(method)
            else:
                # we should check for more domains here when we have them
                pass
        
        # finally, make a domain for the mda - the whole query length
        self.domains.create(
            type='mda',
            full_type_name='Multiple Domain Architecture',
            status_id=0,
            description='Full length query sequence',
            accession='',
            shortname='QUERY',
            ali_from=1,
            ali_to=len(self.fasta_sequence),
            c_evalue=0,
            i_evalue=0,
            homologs_gathered=False,
            intrepid_done=False,
            discern_done=False,
            currently_running_job = 'None'
        )
        return None

    def log_domain_scan_completion(self):
        ''' this function essentially tells the daemon that we are done scanning for the domains '''
        self.domains_scanned = True
        self.save()
        return None

    def pfam_domains(self):
        ''' this function returns the pfam domains found in the query. '''
        return self.domains.filter(type='pfam')

    class Meta:
        db_table = u'protein_analysis_job'

class ProteinDomain(models.Model):
    id = models.AutoField(primary_key=True)
    # what type of domain is this?  pfam, pdb, mda, etc.
    type = models.TextField()
    # holds a full description of this domain type
    full_type_name = models.TextField()
    # foreign key to this job
    protein_analysis_job = models.ForeignKey(ProteinAnalysisJob, related_name='domains')
    # what is the status of this domain?
    status = models.ForeignKey(ProteinDomainStatus, related_name='domain')
    # description of the domain, if available
    description = models.TextField()
    # external accession of the domain
    accession = models.TextField()
    # shortname - for pfam domains, or this could be something else for other domains
    shortname = models.TextField()
    # alignment information
    ali_from = models.IntegerField()
    ali_to = models.IntegerField()
    hmm_from = models.IntegerField()
    hmm_to = models.IntegerField()
    c_evalue = models.FloatField()
    i_evalue = models.FloatField()
    # have we gathered homologs on this domain already?
    homologs_gathered = models.BooleanField()
    # should we run intrepid on this domain?
    intrepid_done = models.BooleanField()
    # had discern been run on this domain?
    discern_done = models.BooleanField()
    # job that is currently running
    currently_running_job = models.TextField()
    
    @property
    def currently_running_job_stages(self):
        ''' This function will report back the percentage done as an integer '''

    class Meta:
        db_table = u'protein_domain'
