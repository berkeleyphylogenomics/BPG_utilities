''' This file contains the api functions for phylobuilder. '''

import re
from piston.handler import BaseHandler
from pfacts003.protein_analysis.models import ProteinAnalysisJob
from pfacts003.protein_analysis.consts import *
from pfacts003.phylofacts.models import OhanaJobs
from datetime import datetime
from django.core.validators import validate_email
from django.core.exceptions import ValidationError
from django.contrib.auth.models import User

class ProteinAnalysis(BaseHandler):
    allowed_methods = ('GET','POST')

    def read(self, request, id=None):
        try:
            job = ProteinAnalysisJob.objects.get(id=id)
        except:
            raise Exception("Job Not Found")
        
        # If the domains have been scanned, the get dictionary will contain a domains category
        if 'domains' in request.GET:
            return dict([( domain.id, {'status': domain.status.status, 'status_bit': domain.status_id, 
                            'currentlyRunningJobStages': domain.currently_running_job_stages, 
                            'type': domain.type} ) for domain in job.domains.all()]) 
                
        else:
            return {
                'id': job.id,
                'status_id': job.status,
                'domains_scanned': job.domains_scanned
            }

    def create(self, request):
        # Create a new Protein Analysis job
        now = datetime.now()
        # This represents today's date.  We use it to find how many jobs the
        # user has submitted today.
        lower_date_bound = datetime(year=now.year, month=now.month, day=now.day)

        # get the ip address from the request
        try:
            ip = request.META['REMOTE_ADDR']
        except:
            ip = '0.0.0.0'

        #######################################################################
        # THROTTLING
        #######################################################################
        
        # is the user logged in?
        if request.user.is_authenticated():
            u = request.user
            # get the jobs submitted by this user today
            user_jobs = OhanaJobs.objects.filter(submitting_user = u, job_type='protein_analysis', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            # should we throttle?
            if len(user_jobs) >= u.profile.get().maximum_phylobuilder_jobs:
                return {
                    'status': 'error',
                    'type': 'throttle',
                    'message': AUTHENTICATED_USER_THROTTLE_MESSAGE % (len(user_jobs))
                }
        else:
            # anonymous user
            # this needs to be changed.  
            u = User.objects.get(username='test')
            # get jobs
            anon_jobs = OhanaJobs.objects.filter(submitting_ip_address = ip, job_type='protein_analysis', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            if len(anon_jobs) >= ANONYMOUS_JOB_LIMIT:
                return {
                    'status': 'error',
                    'type': 'throttle',
                    'message': ANONYMOUS_USER_THROTTLE_MESSAGE % (len(anon_jobs))
                }
        
        ##########################################################################
        # Input processing and FASTA validation
        ##########################################################################
        if ('fasta' in request.POST) and (len(request.POST['fasta']) > 0):
            query = request.POST['fasta']
            # get the header line if it exists
            if query.startswith('>'):
                header = query.split('\n')[0]
            else:
                header = '>QUERY'
            sequence = ''.join(query.split('\n')[1:])
            p = re.compile(r'[\s\n]*')
            fasta_sequence = p.sub('', sequence)
            fasta_sequence = fasta_sequence.upper()
            for (index, aa) in enumerate(fasta_sequence):
                if aa not in ACCEPTED_SEQUENCE_CHARS:
                    return {
                        'status': 'error',
                        'type': 'fasta',
                        'message': 'Unrecognized sequence character %s at position %d' % (aa, index)
                    }
            if len(fasta_sequence) < MINIMUM_INPUT_LENGTH:
                return {
                    'status': 'error',
                    'type': 'fasta',
                    'message': INPUT_SEQUENCE_TOO_SHORT_MESSAGE
                }
            elif len(fasta_sequence) > MAXIMUM_INPUT_LENGTH:
                return {
                    'status': 'error',
                    'type': 'fasta',
                    'message': INPUT_SEQUENCE_TOO_LONG_MESSAGE
                }
        else:
            return {
                'status': 'error',
                'type': 'fasta',
                'message': NO_FASTA_MESSAGE
            }

        ############################################################################################
        # EMAIL VALIDATION
        ############################################################################################
        if ('email' in request.POST) and (len(request.POST['email']) > 0):
            try:
                validate_email(request.POST['email'])
                email = request.POST['email']
            except:
                return {
                    'status': 'error',
                    'type': 'email',
                    'message': EMAIL_INVALID_MESSAGE
                }
        else:
            if request.user.is_authenticated():
                email = request.user.email
            else:
                email = ''

        #############################################################################################
        # JOB COMMENTS VALIDATION
        #############################################################################################
        if ('job-comments' in request.POST):
            if (len(request.POST['job-comments']) == 0):
                job_name = ''
            else:
                if (len(request.POST['job-comments']) <= MAXIMUM_COMMENTS_CHARACTERS):
                    job_name = request.POST['job-comments']
                else:
                    return {
                        'status': 'error',
                        'type': 'job-comments',
                        'message': JOB_COMMENTS_TOO_LONG_MESSAGE
                    }
        else:
            job_name = ''

        ##############################################################################################
        # EMAIL SUBJECT VALIDATION
        ##############################################################################################
        if ('email-subject' in request.POST):
            if len(request.POST['email-subject']):
                if len(request.POST['email-subject']) <= MAXIMUM_EMAIL_SUBJECT_CHARACTERS:
                    email_subject = request.POST['email-subject']
                else:
                    return {
                        'status': 'error',
                        'type': 'email-subject',
                        'message': EMAIL_SUBJECT_TOO_LONG_MESSAGE
                    }
            else:
                email_subject = ''
        else:
            email_subject = ''

        params = {}

        ###############################################################################################
        # PARAMETER VALIDATION
        ###############################################################################################
        # Minimum domain length validation
        if ('minimum-domain-length' in request.POST):
            try:
                minimum_domain_length = int(request.POST['minimum-domain-length'])
            except:
                return {'status':'error', 'type':'domain-length', 'message': MINIMUM_DOMAIN_LENGTH_UNRECOGNIZED_MESSAGE}
            if minimum_domain_length > MAXIMUM_MINIMUM_DOMAIN_LENGTH:
                {'status':'error', 'type':'domain-length', 'message': MINIMUM_DOMAIN_LENGTH_TOO_LONG_MESSAGE}
            elif minimum_domain_length < MINIMUM_MINIMUM_DOMAIN_LENGTH:
                {'status':'error', 'type':'domain-length', 'message': MINIMUM_DOMAIN_LENGTH_TOO_SHORT_MESSAGE}
        else:
            return {'status':'error', 'type':'domain-length', 'message': MINIMUM_DOMAIN_LENGTH_MISSING}

        # domain method validation
        if ('domains-to-scan' in request.POST):
            domain_methods = request.POST['domains-to-scan']

        ###############################################################################################
        # FIX FOR CRAPPY DJANGO PROBLEM
        ###############################################################################################
        try:
            job = ProteinAnalysisJob.objects.get(id=1)
            job.fasta_sequence = 'A'
            job.save()
        except:
            pass
        try:
            job = OhanaJobs.objects.get(id=1)
            job.job_type='protein_analysis'
            job.save()
        except:
            pass

        ###############################################################################################
        # CREATE JOB
        ###############################################################################################
        job = ProteinAnalysisJob.objects.create(
            fasta_header = header,
            fasta_sequence = fasta_sequence,
            #email_subject = email_subject,
            #user_email = email,
            #is_done = False,
            save_job = False,
            minimum_domain_length = minimum_domain_length,
            domain_methods = domain_methods,
            status = 1,
            domains_scanned = False,
            job_path = '',
            #job_name = job_name,
        )

        ###############################################################################################
        # LOG JOB IN OHANA JOB TABLE
        ###############################################################################################
        ojob = OhanaJobs.objects.create(
            submitting_user = u,
            submitting_ip_address = ip,
            job_type = 'protein_analysis'
        )
        job.ohana_job = ojob
        job.save()
        return self.read(None, job.id)

# This class is for submitting a job

class ProteinDomain(BaseHandler):
    allowed_methods = ('GET','POST')
    
    def read(self, request, id):
        ''' This returns the status of the job '''
        return None

    def create(self, request, id):
        ''' Post to this address to continue the analysis '''
               
        return self.read(None, job.id)

