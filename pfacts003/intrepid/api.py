''' This file contains the api functions for intrepid. '''

import re, tempfile, codecs
from piston.handler import BaseHandler
from cStringIO import StringIO
from Bio import AlignIO
from pfacts003.intrepid.models import IntrepidJob, IntrepidStatus
from pfacts003.intrepid.consts import *
from pfacts003.intrepid.utils import get_pdb_file_for_id, get_pdb_fasta_for_id, read_sequence_from_pdb_file, \
                                        mafft, get_chains_from_any_pdb_file, read_sequence_from_any_pdb_file
from pfacts003.phylofacts.models import OhanaJobs
from datetime import datetime
from django.core.validators import validate_email
from django.core.exceptions import ValidationError
from django.contrib.auth.models import User

class Intrepid(BaseHandler):
    allowed_methods = ('GET','POST')

    def read(self, request, id=None):
        try:
            job = IntrepidJob.objects.get(id=id)
        except:
            raise Exception("Job Not Found")

        return job.status_object()

    def create(self, request):
        # Create a new Intrepid job
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
            user_jobs = OhanaJobs.objects.filter(submitting_user = u, job_type='intrepid', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            # should we throttle?
            #if len(user_jobs) >= u.profile.get().maximum_phylobuilder_jobs:
            #    return {
            #        'status': 'error',
            #        'type': 'throttle',
            #        'message': AUTHENTICATED_USER_THROTTLE_MESSAGE % (len(user_jobs))
            #    }
        else:
            # anonymous user
            # this needs to be changed.  
            u = User.objects.get(username='test')
            # get jobs
            anon_jobs = OhanaJobs.objects.filter(submitting_ip_address = ip, job_type='intrepid', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            #if len(anon_jobs) >= ANONYMOUS_JOB_LIMIT:
            #    return {
            #        'status': 'error',
            #        'type': 'throttle',
            #        'message': ANONYMOUS_USER_THROTTLE_MESSAGE % (len(anon_jobs))
            #    }
        ##############################################################################################
        # PDB VALIDATION
        ##############################################################################################
        pdbid = ''
        is_discern = False
        internal_job_name = ''
        pid = ''
        chainid = ''
        chainString = ''
        pdb_file = ''
        has_uploaded_model = False
        if (('pdb-file' in request.POST) and (len(request.POST['pdb-file']) > 0)):
            # start the pdb stuff
            try:
                pdb_file = str(request.POST['pdb-file'])
                has_uploaded_model = True
                if (('pdb' in request.POST) and (len(request.POST['pdb']) > 0)):
                    if ':' in request.POST['pdb']:
                        pid = request.POST['pdb'].split(':')[0].lower()
                        chainid = request.POST['pdb'].split(':')[1].upper()
                    elif '_' in request.POST['pdb']:
                        pid = request.POST['pdb'].split('_')[0].lower()
                        chainid = request.POST['pdb'].split('_')[1].upper()
                    else:
                        pid = request.POST['pdb'][:4].lower()
                        chainid = request.POST['pdb'][4:].upper()
                else:
                    pid = GENERIC_PDB_ID
                if not pdb_file.startswith('HEADER'):
                    pdb_file = ('HEADER    UPLOADED STRUCTURE                       01-JAN-01   %s              \n' % (pid)) + pdb_file
                if not chainid:
                    longest_chain_length = 0
                    this_seq = ''
                    this_tuple_list = []
                    these_chains = get_chains_from_any_pdb_file(pdb_file, pid)
                    for c in these_chains:
                        (se, tl) = read_sequence_from_any_pdb_file(pdb_file, pid, c)
                        if len(se) > longest_chain_length:
                            this_seq = se
                            struct_pos_list = tl
                            chainid = c.strip()
                            chainString = 'No chain supplied. Using chain: %s' % (chainid)
                if (len(this_seq) < MINIMUM_INPUT_LENGTH):
                    return {
                        'status': 'error',
                        'type': 'pdb',
                        'message': 'The longest chain in the structure is %d amino acids.  The minimum for analysis is %s.' % (len(this_seq), MINIMUM_INPUT_LENGTH)
                    }
                fasta_sequence = this_seq
                header = '>lcl|%s|%s Uploaded structure\n' % (pid, chainid)
                internal_job_name = 'Chain %s' % (chainid)
                ali = mafft(header+fasta_sequence+'\n>lcl|FASTA| Sequence from fasta\n'+fasta_sequence)
                fh = StringIO(ali)
                msa = AlignIO.read(fh, 'fasta')
                fh.close()
            except Exception, e:
                print str(e)
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': 'Error parsing file.  Biopython PDB parser cannot recognize file.'
                }    
        elif (('pdb' in request.POST) and len(request.POST['pdb']) > 0):
            is_discern = True
            # Get the file path
            if ':' in request.POST['pdb']:
                pid = request.POST['pdb'].split(':')[0].lower()
                chainid = request.POST['pdb'].split(':')[1].upper()
            elif '_' in request.POST['pdb']:
                pid = request.POST['pdb'].split('_')[0].lower()
                chainid = request.POST['pdb'].split('_')[1].upper()
            else:
                # try to get the chain and pdb id this way
                pid = request.POST['pdb'][:4].lower()
                chainid = request.POST['pdb'][4:].upper()
            pdbid = get_pdb_file_for_id(pid, PDB_LIBRARY_PATH)
            if not pdbid:
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': '%s is not found in the phylofacts pdb structure library' % pid
                }
            (pdb_fasta, internal_job_name, chain_count) = get_pdb_fasta_for_id(pid,chainid,PDB_LIBRARY_PATH)
            if not pdb_fasta:
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': 'Cannot retrieve sequence information for %s' % pid
                }
            if not chainid:
                # store the correct chain we are returning
                chainid = pdb_fasta.split('\n')[0].split()[0].split('|')[2]
                chainString = 'No chain supplied.  Using chain: %s' % chainid
            if chain_count > 1:
                internal_job_name = 'Chains: ' + internal_job_name
            else:
                internal_job_name = 'Chain ' + internal_job_name
            header = pdb_fasta.split('\n')[0]
            fasta_sequence = ''.join(pdb_fasta.split('\n')[1:]) 
            # look for the type of sequence, is it a prot?
            try:
                m = re.search('(?<=mol:)\w+', header)
                pdb_type = m.group(0)
            except:
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': 'Could not ascertain sequence type of %s.' % (pid)
                }                
            if (pdb_type != 'protein'):
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': '%s is of %s type.  Only PDB structures of type protein can be used.' % (pid, pdb_type)
                }
            if (len(fasta_sequence) < MINIMUM_INPUT_LENGTH):
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': '%s is %d amino acids.  The minimum length for analysis is %d' % (pid, len(fasta_sequence), MINIMUM_INPUT_LENGTH)
                }
            # BAD, BAD BAD, BAD BAD BAD BAD, etc.
            try:
                # check the amino acid sequence listed in the pdb file
                (struct_seq, struct_pos_list) = read_sequence_from_pdb_file(pid, chainid, PDB_LIBRARY_PATH)
                # align the sequences from the fasta file and the structure file
                ali = mafft(header+'\n'+fasta_sequence+'\n>lcl|PDB|Sequence from struct\n'+struct_seq)       
                # read the alignment in
                fh = StringIO(ali)
                msa = AlignIO.read(fh, 'fasta')
                fh.close() 
            except:
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': 'Error aligning sequence in structure file with sequence in sequence database'
                }
            # error conditions...don't allow a job creation if any of these are met 
            if ('-' in str(msa[0].seq)) or (len(struct_seq) != len([c for c in str(msa[1].seq) if c.isupper()])):
                return {
                    'status': 'error',
                    'type': 'pdb',
                    'message': 'Error parsing alignment between sequence in structure file and sequence in sequence database.  Possible insertion in sequence from structure file' 
                }
        else:
            return {
                'status': 'error',
                'type': 'pdb',
                'message': 'Either a pdb or an uploaded model is required'
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
 
                    
        ##############################################################################################
        # EMAIL SUBJECT VALIDATION
        ##############################################################################################
        params = {}

        ###############################################################################################
        # PARAMETER VALIDATION
        ###############################################################################################
        homolog_cov = 0.7
        if ('homolog-coverage' in request.POST):
            f = int(request.POST['homolog-coverage'])
            if (f>=0 and f<=99):
                homolog_cov = 1.0*f/100.0
            else:
                return {
                    'status': 'error',
                    'type': 'parameter',
                    'message': 'Minimum column coverage parameter unrecognized'
                }
        seqdb = UNIREF100
        use_old = True
        if ('sequence-database' in request.POST):
            try:
                seqdb = int(request.POST['sequence-database'])
                if seqdb > UNIREF90:
                    use_old = True
                else:
                    use_old = False
                f = SEQ_DBS[seqdb]
            except:
                seqdb = UNIREF100
                use_old = True
        gappiness = 0.5
        if ('minimum-column-coverage' in request.POST):
            try:
                f = float(request.POST['minimum-column-coverage'])
                if (f>=0) and (f<=99):
                    gappiness = 1.0 - (1.0*f/100.0)
                else:
                    return {
                        'status': 'error',
                        'type': 'parameter',
                        'message': 'Minimum column coverage parameter unrecognized'
                    }
            except:
                return {
                    'status': 'error',
                    'type': 'parameter',
                    'message': 'Minimum column coverage parameter unrecognized'
                }
        save_outputs = False
        if ('save-program-outputs' in request.POST):
            try:
                save_outputs = bool(request.POST['save-program-outputs'])
            except:
                save_outputs = False
        maximum_homologs = 2500
        if ('maximum-homologs' in request.POST):
            try:
                maximum_homologs = int(request.POST['maximum-homologs'])
            except:
                maximum_homologs = 2500
        minimum_pfam_dom_len = 40
        if ('minimum-pfam-domain-length' in request.POST):
            try:
                minimum_pfam_dom_len = int(request.POST['minimum-pfam-domain-length'])
                '''if minimum_pfam_dom_len > MAXIMUM_PFAM or minimum_pfam_dom_len < MINIMUM_PFAM:
                    return {
                        'status': 'error',
                        'type': 'parameter',
                        'message': 'Minimum pfam domain length out of bounds'
                    }'''
            except:
                return {
                    'status': 'error',
                    'type': 'parameter',
                    'message': 'Minimun pfam domain length unrecognized'
                }
        hom_iter = 4
        if ('homolog-iterations' in request.POST):
            try:
                hom_iter = int(request.POST['homolog-iterations'])
                '''if hom_iter > MAXIMUM_HOMOLOG_ITERATIONS or hom_iter < MINIMUM_HOMOLOG_ITERATIONS:
                    return {
                        'status': 'error',
                        'type': 'parameter',
                        'message': 'Maximum homolog iterations out of bounds'
                    }'''
            except:
                return {
                    'status': 'error',
                    'type': 'parameter',
                    'message': 'Maximum homolog iterations unrecognized'
                }
        pwid_remote = 0.15
        if ('pwid-query-homolog' in request.POST):
            try:
                f = float(request.POST['pwid-query-homolog'])
                if (f>=0) and (f<=99):
                    pwid_remote = (1.0*f/100.0)
                else:
                    return {
                        'status': 'error',
                        'type': 'parameter',
                        'message': 'Pairwise ID parameter unrecognized'
                    }
            except:
                return {
                    'status': 'error',
                    'type': 'parameter',
                    'message': 'Pairwise ID parameter unrecognized'
                }

        ###############################################################################################
        # FIX FOR CRAPPY DJANGO PROBLEM
        ###############################################################################################
        try:
            job = IntrepidJob.objects.get(id=1)
            job.fasta_sequence = 'A'
            job.save()
        except:
            pass
        try:
            job = OhanaJobs.objects.get(id=1)
            job.job_type='intrepid'
            job.save()
        except:
            pass
        try:
            r = IntrepidJobResidues.objects.get(id=1)
            r.fasta_residue = 'A'
            r.save()
        except:
            pass

        ###############################################################################################
        # CREATE JOB
        ###############################################################################################
        s = IntrepidStatus.objects.get(id=WAITING_FOR_SUBMISSION)        

        job = IntrepidJob.objects.create(
            fasta_header = header,
            fasta_sequence = fasta_sequence,
            user_email = email,
            internal_name = internal_job_name,
            is_done = False,
            has_uploaded_model = has_uploaded_model,
            is_discern_job = True,
            save_job = False,
            status = s,
            max_homologs = maximum_homologs,
            pdb_id = pid+chainid,
            save_program_outputs = save_outputs,
            structure_contains_internal_gaps = False,
            # comment this line out before you commit please.
            #development_job = True,
            development_job = False,
            #job_path = this_directory,
            homolog_evalue = 0.001,
            homolog_coverage_condition = homolog_cov,
            substatus = 'Initialized.\n%s\n' % chainString,
            is_running = False,
            num_iterations = 0,
            num_homologs = 0,
            maximum_divergence = 0.0,
            seq_db = seqdb,
            use_old_blasttools = use_old,
            maximum_homolog_iterations = hom_iter,
            homolog_divergence = pwid_remote,
            gappy_column_condition = gappiness,
            minimum_pfam_length = minimum_pfam_dom_len,
            comparative_model_string = pdb_file
        )

        this_time = datetime.now()
        job.substatus += 'Job %d started at %s\n' % (job.id, this_time.strftime('%a %b %Y, %I:%M:%S %p')) 
        job.substatus += 'Input sequence:\n'
        job.substatus += '%s\n' % (job.fasta_header)
        job.substatus += '%s\n' % (job.pretty_sequence())
        job.substatus += 'Minimum length of Pfam domain for separate analysis: %d\n' % (job.minimum_pfam_length)
        job.substatus += 'ID necessary between query and most remote homolog: %d%%\n' % (int(100 * job.homolog_divergence))
        job.substatus += 'Minimum query (hmm) coverage for inclusion into homolog set: %d%%\n' % (int(100 * job.homolog_coverage_condition))
        job.substatus += 'Maximum iterations of homolog gathering: %d\n' % (job.maximum_homolog_iterations)
        job.substatus += 'Minimum column coverage for inclusion into tree MSA: %d%%\n'% (int(100 * job.gappy_column_condition))
        job.substatus += 'Maximum number of homologs: %d\n' % (job.max_homologs)
        job.substatus += 'Sequence database to gather homologs from: %s\n' % (SEQ_DBS[job.seq_db]['name'])
        job.save()

        ###############################################################################################
        # LOG JOB IN OHANA JOB TABLE
        ###############################################################################################
        ojob = OhanaJobs.objects.create(
            submitting_user = u,
            submitting_ip_address = ip,
            job_type = 'intrepid'
        )
        job.ohana_job = ojob
        job.save()

        if job.is_discern_job:
            # create the intrepid residue objects
            struct_pos = 0
            fali = str(msa[0].seq)
            sali = str(msa[1].seq)
            job.pdb_fasta_alignment_line = fali
            job.pdb_structure_alignment_line = sali
            job.pdb_structure_fasta_offset = struct_pos_list[0][0] - (len(sali) - len(sali.lstrip('-'))) - 1
            job.save()
            for (pos, residue) in enumerate(fali):
                new_res = job.residues.create(
                    fasta_residue = residue,
                    structure_residue = sali[pos],
                    fasta_position = pos + 1,
                    fasta_residue_three_letter = AA_DICT[residue]['3'],
                    discern_rank = 9999,
                    intrepid_rank = 9999,
                    intrepid_cons_js_score = -9999.99,
                    intrepid_global_js_score = -9999.99,
                    discern_score = -9999.99,
                    discern_centrality_score = -9999.99,
                    relative_solvent_accessibility = -9999.99,
                    absolute_solvent_accessibility = -9999.99,
                    global_conservation = -9999.99
                )
                if sali[pos] == '-':
                    new_res.structure_position = 0
                    new_res.structure_residue_three_letter = 'N/A'
                else:
                    new_res.structure_position = struct_pos_list[struct_pos][0]
                    new_res.structure_residue_three_letter = AA_DICT[sali[pos]]['3']
                    struct_pos += 1
                new_res.save()
        return self.read(None, job.id)
