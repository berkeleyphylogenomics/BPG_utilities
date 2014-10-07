# api handlers.py

import sys, hashlib, base64, re, os
from piston.utils import throttle
from piston.handler import BaseHandler
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Seq import Seq
from bpg.Classifiers import HMMBLASTSequenceClassifier
from pfacts003.phylofacts.models import UniProt, UniProtIDMapping, Pfam, Family, FatcatJob, FatcatJobFamily, FxnSitePredictionJob, FxnSitePredictionJobStatus, MetacycReaction, TreeNode, PhyloFactsUserProfile, OhanaJobs, UniProtGene
from pfacts003.fatcat.models import FatcatJob as Fatcat2Job
from pfacts003.fatcat.models import FatcatJobFamily as Fatcat2JobFamily
from pfacts003.fatcat.consts import *
from pfacts003.utils.extract_sequence_from_fasta import extract_sequence_from_fasta
from pfacts003.utils.phylofacts_email import send_email_from_phylofacts
from pfacts003.phylofacts.forms import NewSequenceSearchForm
import operator
from django.db.models import Q
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.models import User
from pysolr import Solr
from pfacts003.utils.id_patterns import fasta_re, is_uniprot_identifier_format, is_uniprot_accession_format, bpgid_re, phog_re
from django.core.validators import validate_email
from django.core.exceptions import ValidationError
from django.db.utils import IntegrityError
from django.utils.cache import add_never_cache_headers
from datetime import datetime, timedelta
import pprint
try:
   import json #json isn't in 2.4
except ImportError:
   import simplejson as json

class PhyloFactsPasswordChangeHandler(BaseHandler):
    allowed_methods = ('POST')

    def create(self, request):
        if request.user.is_authenticated():
            try:
                u = request.user
                u.set_password(request.POST['password'])
                u.save()
                return {'success':'true'}
            except:
                return {'error':'true'}
        else:
            return {'error':'true'}

class PhyloFactsPasswordResetHandler(BaseHandler):
    allowed_methods = ('POST')

    def create(self, request):
        email_to = request.POST['email']
        try:
            u = User.objects.get(email=email_to)
        except:
            return {'error':'email', 'message': 'Could not find unique user with this email.  <a href="/phylofacts/user_create/">Click here</a> to create a new user.'}
        new_password = User.objects.make_random_password()
        try:
            u.set_password(new_password)
        except:
            return {'error': 'set', 'message': 'Error setting new password'}
        send_email_from_phylofacts(email_to, "This is an automatically generated password to your PhyloFacts user account.  Use it to log in to phylofacts at: http://makana.berkeley.edu/phylofacts/login/.\n\nUsername: %s\nPassword: %s" % (u.username, new_password), "This is an automatically generated password to your PhyloFacts user account.  Use it to log in to <a href='http://makana.berkeley.edu/phylofacts/login/'>PhyloFacts</a><br /><br />Username: %s<br />Password: %s" % (u.username, new_password), "Your PhyloFacts password has been reset.")
        return {'success': '', 'message': 'An email has been sent to %s containing your new password' % email_to}

class PhyloFactsUserCreationHandler(BaseHandler):
    allowed_methods = ('POST')

    #@throttle(2, 24*60*60)
    def create(self, request):
        try:
            validate_email(str(request.POST['email'])[:254])
            email = str(request.POST['email'])[:254]
        except ValidationError:
            return {'error': 'email', 'message': 'Please input a valid email address.'}
        if (len(request.POST['first_name']) > 0):
            first_name = request.POST['first_name'][:30]
        else:
            return {'error': 'first_name','message': 'First name required.'}
        if (len(request.POST['last_name']) > 0):
            last_name = request.POST['last_name'][:30]
        else:
            return {'error':'last_name', 'message':'Last name required.'}
        if (len(request.POST['username']) == 0):
            return {'error':'username', 'message':'Username is required.'}
        if (len(request.POST['username']) > 30):
            return {'error':'username', 'message':'Username must be fewer than 30 characters.'}
        institution = request.POST['institution_name'][:50]
        institution_country = request.POST['institution_country']
        institution_address = str(request.POST['institution_address'])[:200]
        if len(request.POST['institution_type']) < 2:
            institution_type = request.POST['institution_other'][:30]
        else:
            institution_type = request.POST['institution_type']
        if len(request.POST['user_position']) < 2:
            user_position = request.POST['user_position_other'][:50]
        else:
            user_position = request.POST['user_position']
        # create user
        try:
            u = User.objects.get(username='test')
            v = PhyloFactsUserProfile.objects.get(user = u)
            v.institution = 'test' 
            u.is_staff = False
            password = User.objects.make_random_password()
            username = request.POST['username']
            User.objects.create_user(username, email=email, password=password)
            u = User.objects.get(username=username)
            u.first_name = first_name
            u.last_name = last_name
            u.save()
            u.profile.create(institution = institution, institution_address = institution_address, institution_country = institution_country, user_position = user_position, maximum_fatcat_jobs = 20)
            send_email_from_phylofacts(email, "This is an automatically generated password to your PhyloFacts user account.  Use it to log in to phylofacts at: http://makana.berkeley.edu/phylofacts/login/.\n\nUsername: %s\nPassword: %s" % (u.username, password), "This is an automatically generated password to your PhyloFacts user account.  Use it to log in to <a href='http://makana.berkeley.edu/phylofacts/login/'>PhyloFacts</a><br /><br />Username: %s<br />Password: %s" % (u.username, password), "Your PhyloFacts user account has been created.") 
        except IntegrityError:
            return {'error':'email','message':'This username already exists.  Forgot your password?  <a href="/phylofacts/password_reset">Click here</a>'}
        except:
            return {'error':'user_create','message':'There was an error in user creation, please notify the PhyloFacts webmaster by email at phylofacts.webmaster@gmail.com'}
        return {'message': 'User created.  A temporary password was emailed to %s.' % email}

class PhyloFactsLogoutHandler(BaseHandler):
    allowed_methods = ('POST')

    def create(self, request):
        logout(request)
        return {'logout':'yes'}    

class PhyloFactsLoginHandler(BaseHandler):
    allowed_methods = ('POST')

    def create(self, request):
        try:
            username = request.POST['username']
            password = request.POST['password']
        except:
            return {'error': True}
        # logout if we are already logged in
        if request.user.is_authenticated():
            logout(request)
        # try to authenticate the user
        user = authenticate(username=username, password=password)
        if user is not None:
            if user.is_active:
                p = user.profile.get()
                if not p.has_ever_logged_in:
                    p.has_ever_logged_in = True
                    p.save() 
                login(request, user)
                return {'success': True}
            else:
                return {'error': True}
        else:
            return {'error': True}
        
class PhyloScopeAJAXRequestHandler(BaseHandler):
    allowed_methods = ('POST')

    def create(self, request):
        if ("tree_node_id" in request.POST):
            try:
                node = TreeNode.objects.get(id = int(request.POST["tree_node_id"]))
            except:
                return {'status' : 1, 'message' : 'No matching tree node id'}
            try:
                d = node.consensus_description.get().consensus_uniprot_description
                t = node.mrca.get().mrca
            except:
                return {'status' : 1, 'message' : 'Could not find tree node information'}
            ret_obj = {}
            if d:
                ret_obj['consensus_description'] = d
            if t:
                ret_obj['scientific_name'] = t.scientific_name
                ret_obj['taxon_id'] = t.id
                if t.common_name:
                    ret_obj['common_name'] = t.common_name
            ret_obj['status'] = 0
            return ret_obj
        else:
            return {'status': 1, 'message': 'No field tree_node_id in POST'}

class PhyloFactsSearchHandler(BaseHandler):
    # This function handles the search box on our phylofacts website.  
    # There is a string "input" field in the post to this handler, this
    # function checks if it is in various formats, before it flags an error
    # and returns that.  The 'status' field in the return must be there, and
    # the parser on the webform handles the accompanying fields according to that
    # field.  Now, the only 2 status options parsed there are redirect or error.
    # To preserve the separate error messages, I added a searchtype input field.
    # This allows the error messages to be sent to the right div.
    allowed_methods = ('POST')

    def create(self, request):
        if (("input" in request.POST) and ("div_id" in request.POST)):
            query = request.POST["input"]
            displaytodiv = request.POST["div_id"]
            """General Accession and Identifier Search
                Right now we support:
                    bpg0240116
                    APAF_HUMAN
                    P30559
                    PHOG0093626_01452
                    PF12345
                    Pfam short name (Pfam.name)
             """
            if bpgid_re.match(query.lower()):
                try:
                    fam = Family.find_by_family_id(int(query.lower().strip("bpg")))
                except Family.DoesNotExist:
                    return {'status' : 'error', 'displaytodiv' : displaytodiv, 
                            'message' : '%s is not a valid bpg accession' % query} 
                return {'status' : 'redirect', 'url': '/phylofacts/family/%s' % query.lower()}
            # TODO do the phog matching stuff
            if phog_re.match(query):
                return {'status' : 'redirect', 'url': '/phog/%s' % query}

            if is_uniprot_identifier_format(query) or is_uniprot_accession_format(query):
                try:
                    uni = UniProt.find_by_accession_or_identifier(query.upper())
                except UniProt.DoesNotExist:
                    return {'status' : 'error', 'displaytodiv' : displaytodiv,
                            'message' : '%s is not included in the PhyloFacts database' % query.upper()}    
                return {'status' : 'redirect', 'url': '/phylofacts/sequence/UniProt/%s' % query.upper()}

            pfam_acc = re.match("PF\d{5}", query.upper())
            if pfam_acc is not None:
                q = Pfam.objects.filter(accession=pfam_acc.group(0))
                if len(q) != 0:
                    return {'status' : 'redirect', 'url': '/phylofacts/pfam/%s' % pfam_acc.group(0)}
                else:
                    return {'status' : 'error', 'displaytodiv' : displaytodiv,
                            'message' : '%s is not a PhyloFacts-Pfam accession' % pfam_acc.group(0)}
            q = Pfam.objects.filter(name__iexact = query.lower())
            if len(q) != 0:
                return {'status' : 'redirect', 'url': '/phylofacts/pfam/%s' % q[0].accession}


            # Both the rxn_acc below are for finding biocyc reactions. TODO: Combine this and other queries.
            rxn_acc = re.match("\d+\.[\d+|-]\.[\d+|-]\.[\d+|-]", query.upper())
            if rxn_acc is not None:
                q = MetacycReaction.objects.filter(ec_number=rxn_acc.group(0))
                if len(q) != 0:
                    return {'status' : 'redirect', 'url': '/phylofacts/biocyc_rxn/%s' % q[0].id}
                else:
                    return {'status' : 'error', 'displaytodiv' : displaytodiv,
                            'message' : '%s is not a PhyloFacts-Pfam accession' % q[0].id}
            rxn_acc = re.match("[\w\-]{0,100}RXN[\w\-]{0,100}", query.upper())
            if rxn_acc is not None:
                q = MetacycReaction.objects.filter(id=rxn_acc.group(0))
                if len(q) != 0:
                    return {'status' : 'redirect', 'url': '/phylofacts/biocyc_rxn/%s' % q[0].id}
                else:
                    return {'status' : 'error', 'displaytodiv' : displaytodiv,
                            'message' : '%s is not a PhyloFacts-Pfam accession' % q[0].id}
            # Attempt to look up by other IDs (kind of dumb right now, it just takes the first hit)
            uniprot = UniProtIDMapping.find_uniprot(query)
            if uniprot:
                return {'status' : 'redirect', 'url': '/phylofacts/sequence/UniProt/%s' % uniprot.accession}


         # nothing matched, return a redirect to go to the solrsearch
            return {'status' : 'redirect', 'url' : '/phylofacts/search/?query=%s' % query}
        else:
            return {'status' : 'error', 'displaytodiv' : '#general-errors', 'message' : 'No input.'}

class HMMBLASTSequenceClassifierHandler(BaseHandler):
    allowed_methods = ('GET', 'POST')
    
    def read(self, request, id):
        # Pretty much just returns a status request.
        return HMMBLASTSequenceClassifier.get_job_by_id(id).status()

    def create(self, request):
        force = False
        form = NewSequenceSearchForm(request.POST)
        if "f" in request.POST:
            force = (request.POST["f"] == "T")
        if form.is_valid():
            input = form.cleaned_data['input']
            defline, sequence, errors = extract_sequence_from_fasta(input)
            if defline == '' or defline is None:
                defline = ">Sequence without header submitted."
        else:
            return {'status' : 'error', 'errors' : form.errors['input']}
        
        uniprot_object = []
        if not force:
            uniprot_object = UniProt.find_by_sequence(sequence)
        if uniprot_object:
            return {'status' : 'existing', 'accession' : str(uniprot_object.accession)}
        else:
            classifierjob = HMMBLASTSequenceClassifier.create_job(defline, sequence)
            classifierjob.run_HMM_BLAST_classification()
            return {'status' : 'created', 'id' : classifierjob.get_job_id()}


'''Provides a rest handler to a fairly abritrary sequence search

A JSON blob of the search parameters is passed in, and we return a list of uniprot accessions

The interface is inspired by PlasmoDB (or will be anyways)

The basic overview is that we are building a limited ORM with a REST API.  We'll probably do this stupidly, but that's okay

I tried doing this simply, copying the quite wonderful DBIx::Class, but then I ran into the restrictions of both python and javascript,
both of whom strongly fought against all my stupid ideas.

The most inside object is a filter clause:

    {
        field: "accession",
        operator: "exact",
        value: "B4KB41"
    }

These are wrapped inside "Sets"

    {
        operator: "AND",
        children: [
            ...
            ...
        ]
    }

'''

# THIS TOTALLY NEEDS TO BE SPLIT INTO A MODULE
class SequenceQueryClause:
    def __init__(self, obj):
        self.field = obj["field"]
        self.operator = obj["operator"]
        self.value = obj["value"]

        # This isn't pretty.  Sorry
        if self.field == 'go_exp_evidence':
            self.field = 'go_annotations__go_evidence__priority'
            if self.value == "0":
                self.operator = 'not_lte'
            else:
                self.operator = 'lte'
            self.value = 5

    def toQ(self):
        if self.operator.startswith( 'not_' ):
            matches = re.search('not_(.*)', self.operator)
            return ~Q( ( self.field + "__" + matches.group(1), self.value  ))
        else: 
            return Q( ( self.field + "__" + self.operator, self.value  ))

class SequenceQuerySet:
    def __init__(self, obj):
        self.operator = obj["operator"]
        self.children = []

        for child in obj["children"]:
            if "value" in child:
                self.children.append( SequenceQueryClause(child) )
            else:
                self.children.append( SequenceQuerySet(child) )


    def  toQ(self):
        if self.operator == "OR":
            op = operator.or_
        else:
            op = operator.and_
        return reduce(op, map(lambda c: c.toQ(), self.children) )


class SequenceQuery:
    '''This is a wrapper for the previous two classes.

       You give it a filter dict object thingy, and it gives you back a set of UniProts

       well, it more shoves them into self.query_set, because fetching them all from the DB
       would just be silly.

       What sucks about all this is that we go from filtering, to collecting, but we should be able
       to do this all this in any point of nesting.  It'd be nice if IN clauses were fast, but they force
       join order.


       For right now, maybe we limit collect orthologs to a top level thing.

       That makes me feel kind of dirty, but I know how to do it, so I guess it works for now.

       I guess for right now I'll assume that every filter has an orthologs() between it.  This is probably limiting and stupid
       So, if you are here two years from now, and you think I'm an idiot for this, I owe you a beer.
    '''

    def __init__(self, filters):
        self.Qs = []
        for filter in filters:
            self.Qs.append( SequenceQuerySet(filter).toQ() )

        self.query_set = UniProt.objects

        first = True
        for Q in self.Qs:
            if not first:
                self.query_set = self.query_set.orthologs(0) 
            self.query_set = self.query_set.filter(Q)
            first = False

        # This is really cracked out, I'm sorry.  The way django's orm works is a bit.. off
        if len(filters) > 1:
            self.query_set = self.query_set.extra(select={"phog_accession": "phog_membership.phog_accession"})
            self.query_set = self.query_set.extra(select={"phog_uniprot_accession": 'T5.accession'})
            self.query_set = self.query_set.extra(select={"phog_uniprot_de": 'T5.de'})


class SequenceHandler(BaseHandler):
    allowed_methods = ('GET')

    def __init__(self):
        self.filters = []

    def read(self, request): # this should take an id to do something smart later.  not now.  -- I miss chained methods!!
        filters = json.loads(request.GET['filters'])
        total = 14295764 #UniProt.objects.count() - it's silly to do this each time - it will be fast in 9.2 (but still silly)


        if "iDisplayStart" in request.GET:
            start = int(request.GET["iDisplayStart"])
        else:
            start = 0

        if "iDisplayLength" in request.GET:
            length = int(request.GET["iDisplayLength"])
        else:
            length = 10

        if "sEcho" in request.GET:
            sEcho = request.GET["sEcho"]
        else:
            sEcho = "1"

        # oh you silly scraper you
        if length > 100:
            length = 100

        end = start + length
        try:
            sequence_query = SequenceQuery(filters)
        except:
            # instead of a 404, we bail here if we don't have a filter
            return {
                "iTotalRecords": total,
                "iTotalDisplayRecords": 0,
                "aaData": [],
                "sEcho": sEcho,
            }

        #count = sequence_query.query_set.count()
        count = 0


        # TODO Respect the params
        if len(filters) > 1:
            fields = ["accession", "uniprot_identifier", "de", "phog_accession", "phog_uniprot_accession", "phog_uniprot_de" ] 
        else:
            fields = ["accession", "uniprot_identifier", "de"] 
        data = sequence_query.query_set.order_by("accession").values(*fields)[start:end]

        
        return {
            "iTotalRecords": total,
            "iTotalDisplayRecords": count,
            "aaData": data,
            "sEcho": sEcho,
            "SQL" : sequence_query.query_set.query
        }


# Wrapper class to convert our query form stuff into solr requests, and then 
# shuffle it between solr and data tables
class SolrHandler(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request):
        if "sSearch" in request.GET:
            query = request.GET["sSearch"]
        else:
            query = "*:*"

        if "iDisplayStart" in request.GET:
            start = int(request.GET["iDisplayStart"])
        else:
            start = 0

        if "iDisplayLength" in request.GET:
            length = int(request.GET["iDisplayLength"])
        else:
            length = 10

        if "sEcho" in request.GET:
            sEcho = request.GET["sEcho"]
        else:
            sEcho = "1"

        # oh you silly scraper you
        if length > 100:
            length = 100

        results = self.conn.search(
            query,
            start=start,
            rows=length
            )

        if "sEcho" in request.GET:
            sEcho = request.GET["sEcho"]
        else:
            sEcho = "1"

        return {
            "iTotalRecords": self.total_records,
            "iTotalDisplayRecords": results.hits,
            "aaData": results.docs,
            "sEcho": sEcho,
        }

class SolrUniProtHandler(SolrHandler):
    def __init__(self):
        self.conn = Solr('http://10.0.0.2:8983/solr/uniprot/')
        self.total_records = 14295764

class SolrFamilyHandler(SolrHandler):
    def __init__(self):
        self.conn = Solr('http://10.0.0.2:8983/solr/family/')
        self.total_records = 92800 

class SolrPhogHandler(SolrHandler):
    def __init__(self):
        self.conn = Solr('http://10.0.0.2:8983/solr/phog/')
        self.total_records = 2643116

class FunctionalSitePredictor(BaseHandler):
    allowed_methods = ('GET', 'POST')

    def read(self, request, jobID):
        job = FxnSitePredictionJob.objects.get(id=jobID)
        if not job:
            raise Exception("Job Not Found")
        
        return {
            'id': job.id,
            'status_id': job.status_id,
            'status': job.status.status,
            'num_homologs':job.num_homologs
        }

    def create(self, request):
        if "fasta" in request.POST:
            #try:
            raw = request.POST["fasta"].split('\n')
            header = raw.pop(0)
            body = ''.join(raw)

            if header[0] != ">":
                body = header+''.join(raw)
                header = ">lcl|PF_FSP_JOB_Q|Sequence with no header submitted."

            header = header[1:] # strip the >

            seq = Seq(body, IUPACProtein)
            try:
                num_iterations = int(request.POST["jackhmmer_iterations"])
                pid = (float(request.POST["treecut_pid"])*0.01)
                evalue = float(request.POST["jackhmmer_evalue"])
                email = request.POST["email"]
            except:
                raise
            # TODO Add in validation.  Because BioPython Sucks

            newjob = FxnSitePredictionJob.objects.create(
                    fasta_header	    = header,
                    fasta_sequence      = str(seq),
                    jackhmmer_iterations= num_iterations,
                    treecut_pid         = pid,
                    jackhmmer_evalue    = evalue,
                    user_email          = email,
                    status_id	        = 1
            )
            newjob_id = FxnSitePredictionJob.objects.all().order_by('id').reverse()[0].id
            return self.read(None, newjob_id)
            #except:
            #    raise Exception("Input is incorrect, check FASTA and parameters.")

# This is the master part of the fatcat jobs
# External users can POST new jobs, and GET
# previous ones
class FatCat(BaseHandler):
    allowed_methods = ('GET', 'POST')

    def read(self, request, fatcat_id=None):
        job = FatcatJob.objects.get(id=fatcat_id)
        if not job:
            raise Exception("Job Not Found")

        return {
            'id': job.id,
            'status_id': job.status_id,
            'status': job.status.status,
            'length': len(job.fasta_sequence),
        }

    def create(self, request):
        now = datetime.now()
        lower_date_bound = datetime(year=now.year, month=now.month, day=now.day)
        try:
            ip = request.META['REMOTE_ADDR']
        except:
            ip = "0.0.0.0"
        if request.user.is_authenticated():
            u = request.user
            # throttle for anonymous user
            user_jobs = OhanaJobs.objects.filter(submitting_user = u, job_type='fatcat', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            if len(user_jobs) >= u.profile.get().maximum_fatcat_jobs:
                return {
                    'status': 'error', 
                    'type': 'throttle', 
                    'message': 'You have submitted %d jobs since midnight PST.  You must wait until tomorrow before submitting again.  Please email phylofacts.webmaster@gmail.com to increase your daily job submission quota.' % (len(user_jobs))
                    }            
        else:
            u = User.objects.get(username='test')
            # throttle for anonymous user
            anon_jobs = OhanaJobs.objects.filter(submitting_ip_address = ip, job_type='fatcat', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            if len(anon_jobs) >= 20:
                return {
                    'status': 'error',
                    'type': 'throttle',
                    'message': 'We have received %d FAT-CAT submissions from this IP address in the last 24 hours.  This is the maximum allowed for non-registered users.  To become a registered user and increase your quota, please go <a href="/phylofacts/user_create/">here</a>.  Otherwise, you can submit new jobs after midnight PST.' % (len(anon_jobs))
                    }
        if ('header' in request.POST):
            if len(request.POST['header']) > 150:
                return {'status':'error', 'type':'header', 'message':'Header must be &lt;150 characters'}
            if len(request.POST['header']) == 0:
                seq_header = ''
            else:
                seq_header = request.POST['header']
        else:
            seq_header = ''
        if ("fasta" in request.POST) and (len(request.POST['fasta']) > 0):
            p = re.compile(r'[\s]*')
            fasta = p.sub('', request.POST['fasta'])
            if fasta[:1] == '>':
                return {'status': 'error', 'type': 'fasta', 'message': '<span class="tip" title="Example: MTNIRKSHPLMKIVNDAFIDLPAPSNISW">FASTA sequence input not accepted.<br />Please input raw sequence.</span>'}
            fasta = fasta.upper()
            achars = 'ABCDEFGHIJKLMNOPQRSTUVWYXZ'
            nondna = 'BDEFHIJKLMNOPQRSVWYXZ' # no GCATU
            for (index, aa) in enumerate(fasta):
                if aa not in achars:
                    return {
                        'status': 'error',
                        'type': 'fasta',
                        'message': 'Error %s at character %d' % (index)
                           }
            is_dna = True
            for aa in fasta:
                if aa in nondna:
                    is_dna = False
                    break
            if is_dna:
                return {
                    'status': 'error',
                    'type': 'fasta',
                    'message': 'Sequence appears to be nucleotide.<br />Only protein sequences are supported.'
                       }
            if (len(fasta) < 15):
                return {'status': 'error', 'type': 'fasta', 'message': 'Sequence too short.<br />FAT-CAT does not accept inputs<br />&lt;15 amino acids'}
            if (len(fasta) > 2000):
                return {'status': 'error', 'type': 'fasta', 'message': 'Sequence too long.<br />FAT-CAT does not accept inputs<br />&gt; 2000 amino acids'}
 
        else:
            return {'status':'error', 'type': 'fasta', 'message': 'FASTA input required'}
        if ("email" in request.POST) and (len(request.POST['email']) > 0):
            try:
                validate_email(request.POST['email'])
                email = request.POST['email']
            except:
                return {'status':'error', 'type':'email', 'message': 'Invalid email'}
        else:
            if request.user.is_authenticated():
                email = request.user.email
            else:
                email = ""
        if ("job-name" in request.POST):
            if (len(request.POST['job-name']) == 0):
                job_name = ''
            else:
                if (len(request.POST['job-name']) < 200):
                    job_name = request.POST['job-name']
                else:
                    return {'status':'error', 'type':'job-name', 'message': 'Comments must be &lt;200 characters'}
        else:
            job_name = ""
        if ('email-subject' in request.POST):
            if (len(request.POST['email-subject']) == 0):
                email_subject = ''
            else:
                if (len(request.POST['email-subject']) < 100):
                    email_subject = request.POST['email-subject']
                else:
                    return {'status':'error', 'type':'email-subject', 'message': 'No more than 100 characters'}

        try:
            cgfev = float(request.POST['stage-one-eval'])
            cgfmqc = float(request.POST['stage-one-mda-qcov'])/100
            cgfmhc = float(request.POST['stage-one-mda-hcov'])/100
            cgfphc = float(request.POST['stage-one-pfam-hcov'])/100
            cgbnev = float(request.POST['stage-two-eval'])
            cgbnmqc = float(request.POST['stage-two-mda-qcov'])/100
            cgbnmhc = float(request.POST['stage-two-mda-hcov'])/100
            cgbnphc = float(request.POST['stage-two-pfam-hcov'])/100
            cgbnqtshp = float(request.POST['query-hmm-pid'])
            cgbnsmp = float(request.POST['subtree-min-pid'])
            cecom = request.POST['ec-checkboxes']
            ceckt = int(request.POST['ec-kerf-threshold'])
            cecomr = int(request.POST['ec-orthology-number'])
            ccs = int(request.POST['cluster-similarity'])
            coc = int(request.POST['ortholog-coverage'])
            cqc = int(request.POST['query-coverage'])
            cmpfo = int(request.POST['minimum-pid-for-orthology'])
            if (request.user.is_authenticated() and request.user.is_staff):
                ccudbp = float(request.POST['consensus-uniprot-description-base-parameter'])/100
                ccutfhc = float(request.POST['consensus-uniprot-threshold-high'])/100
                ccutfmc = float(request.POST['consensus-uniprot-threshold-medium'])/100
            else:
                ccudbp = 0.66
                ccutfhc = 0.7
                ccutfmc = 0.4
        except:
            return {'status': 'error', 'type': 'parameter', 'message': 'Malformed parameters.  Please check program parameters below and submit again.'}
        # print u.username
        # This is a terrible fix for a documented problem with django.  We should upgrade django.
        # Until then, do this.
        # see https://groups.google.com/forum/?fromgroups=#!topic/django-developers/YbOHeu3s7vo
        try:
            job = FatcatJob.objects.get(id=1000)
            job.fasta = 'A'
            job.save()
        except:
            pass
        try:
            oj = OhanaJobs.objects.get(id=1)
            oj.job_type='fatcat'
            oj.save()
        except:
            pass
        job = FatcatJob.objects.create(
                fasta_sequence = fasta,
                status_id = 1,
                user_email = email,
                criteria_get_families_e_value = cgfev,
                criteria_get_families_mda_query_coverage = cgfmqc,
                criteria_get_families_mda_hmm_coverage = cgfmhc,
                criteria_get_families_pfam_hmm_coverage = cgfphc,                
                criteria_get_best_nodes_e_value = cgbnev,
                criteria_get_best_nodes_mda_query_coverage = cgbnmqc,
                criteria_get_best_nodes_mda_hmm_coverage = cgbnmhc,
                criteria_get_best_nodes_pfam_hmm_coverage = cgbnphc,                
                criteria_get_best_nodes_query_to_subtree_hmm_pid = cgbnqtshp,
                criteria_get_best_nodes_subtree_min_pid = cgbnsmp,
                criteria_enclosing_clade_orthology_methods = cecom,
                criteria_enclosing_clade_kerf_threshold = ceckt,
                criteria_enclosing_clade_orthology_methods_required = cecomr,
                criteria_cluster_similarity = ccs,
                criteria_ortholog_coverage = coc,
                criteria_query_coverage = cqc,
                criteria_consensus_uniprot_base_parameter = ccudbp,
                criteria_consensus_uniprot_threshold_for_high_confidence = ccutfhc,
                criteria_consensus_uniprot_threshold_for_medium_confidence = ccutfmc,
                criteria_minimum_pid_to_query_for_orthology = cmpfo,
                submitter_ip_address = ip,
                submitter_user_account = u,
                fatcat_job_name = job_name,
                email_subject = email_subject
        )
        if seq_header:
            header = seq_header
        else:
            header = "FAT-CAT job %d" % job.id
        job.fasta_header = header
        job.save()
        # log this job in the job table.
        ojob = OhanaJobs.objects.create(
            submitting_user = u,
            submitting_ip_address = ip,
            job_type = 'fatcat'
        )
        # reference the ohana job entry in the fatcat job table
        job.ohana_job = ojob
        job.save()
        return self.read(None, job.id)
        
        
# A lot of this belongs in the model.
from pfacts003.utils.annotation import TreeAnnotation, GOAnnotation, summarize, consensus
class FatCatSummary(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        job = FatcatJob.objects.get(id=fatcat_id)
        if not job:
            raise Exception("Job Not Found")

        if job.status_id < 9:
            return "Not Yet Available"

        return job.summary()

class FatCatResults(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        job = FatcatJob.objects.get(id=fatcat_id)
        if not job:
            raise Exception("Job Not Found")
        if job.status_id < 9:
            return "Not Yet Available"
        return job.results()

class FatCatOrthologs(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        job = FatcatJob.objects.get(id=fatcat_id)
        if not job:
            raise Exception("Job Not Found")
        if job.status_id < 9:
            return "Not Yet Available"
        return [o.accession for o in job.orthologs() ]

# New FAT-CAT begins here

class Fatcat2(BaseHandler):
    allowed_methods = ('GET','POST')

    def read(self, request, fatcat_id):
        job = Fatcat2Job.objects.get(id=fatcat_id)
        if not job:
            raise Exception('Job Not Found')

        if not job.is_done:
            now = datetime.now()
            td = now - job.created_at
            # this is weird...
            total_seconds_int = int((float(td.microseconds + (td.seconds + td.days * 24 *3600) * 10**6)) / 10**6)
            if job.status_id >= 7:
                return {
                    'id': job.id,
                    'job-run-time': '%d seconds' % total_seconds_int,
                    'stage1-families': len(job.get_families_passing_stage1()),
                    'stage2-families': len(job.get_families_passing_stage2()),
                    'status_id': job.status_id,
                    'status': job.status.status
                }
            elif job.status_id >= 5:                
                return {
                    'id': job.id,
                    'job-run-time': '%d seconds' % total_seconds_int,
                    'stage1-families': len(job.get_families_passing_stage1()),
                    'status_id': job.status_id,
                    'status': job.status.status
                }
            return {
                'id': job.id,
                'job-run-time': '%d seconds' % total_seconds_int,
                'status_id': job.status_id,
                'status': job.status.status
            }
        else:
            return {
                'id':job.id,
                'status_id': job.status_id,
                'job-run-time': '',
                'status': job.status.status,
                'done': True
                }

    def create(self, request):
        # Create a new FAT-CAT job here
        # implement some kind of job throttling
        now = datetime.now()
        # Today's date.
        lower_date_bound = datetime(year=now.year, month=now.month, day=now.day)
        # get the ip address from the request
        try:
            ip = request.META['REMOTE_ADDR']
        except:
            ip = '0.0.0.0'
        
        if request.user.is_authenticated():
            u = request.user
            # throttle for logged in user
            user_jobs = OhanaJobs.objects.filter(submitting_user = u, job_type='fatcat2', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            if len(user_jobs) >= u.profile.get().maximum_fatcat_jobs:
                return {
                    'status': 'error',
                    'type': 'throttle', 
                    'message': 'You have submitted %d jobs since midnight PST.  You must wait until tomorrow before submitting again.  Please email phylofacts.webmaster@gmail.com to increase your daily job submission quota.' % (len(user_jobs))
                    }
        else:
            u = User.objects.get(username='test')
            # throttle for anonymous user
            anon_jobs = OhanaJobs.objects.filter(submitting_ip_address = ip, job_type='fatcat2', job_create_time__range = [lower_date_bound, now]).order_by('job_create_time')
            if len(anon_jobs) >= 20:
                return {
                    'status': 'error',
                    'type': 'throttle',
                    'message': 'We have received %d FAT-CAT submissions from this IP address in the last 24 hours.  This is the maximum allowed for non-registered users.  To become a registered user and increase your quota, please go <a href="/phylofacts/user_create/">here</a>.  Otherwise, you can submit new jobs after midnight PST.' % (len(anon_jobs))
                    }
        if ('header' in request.POST):
            if len(request.POST['header']) > 150:
                return {'status':'error', 'type':'header', 'message':'Header must be &lt;150 characters'}
            if len(request.POST['header'].strip().lstrip(">").strip()) == 0:
                seq_header = ''
            else:
                seq_header = request.POST['header'].strip().lstrip(">").strip()
        else:
            seq_header = ''
        if ("fasta" in request.POST) and (len(request.POST['fasta']) > 0):
            p = re.compile(r'[\s]*')
            fasta = p.sub('', request.POST['fasta'])
            if fasta[:1] == '>':
                return {'status': 'error', 'type': 'fasta', 'message': '<span class="tip" title="Example: MTNIRKSHPLMKIVNDAFIDLPAPSNISW">FASTA sequence input not accepted.<br />Please input raw sequence.</span>'}
            fasta = fasta.upper()
            achars = 'ABCDEFGHIJKLMNOPQRSTUVWYXZ'
            nondna = 'BDEFHIJKLMNOPQRSVWYXZ' # no GCATU
            for (index, aa) in enumerate(fasta):
                if aa not in achars:
                    return {
                        'status': 'error',
                        'type': 'fasta',
                        'message': 'Error %s at character %d' % (index)
                           }
            is_dna = True
            for aa in fasta:
                if aa in nondna:
                    is_dna = False
                    break
            if is_dna:
                return {
                    'status': 'error',
                    'type': 'fasta',
                    'message': 'Sequence appears to be nucleotide.<br />Only protein sequences are supported.'
                       }
            if (len(fasta) < 15):
                return {'status': 'error', 'type': 'fasta', 'message': 'Sequence too short.  FAT-CAT does not accept inputs &lt;15 amino acids'}
            if (len(fasta) > 2000):
                return {'status': 'error', 'type': 'fasta', 'message': 'Sequence too long.  FAT-CAT does not accept inputs &gt; 2000 amino acids'}
        else:
            return {'status':'error', 'type': 'fasta', 'message': 'FASTA input required'}
        if ("email" in request.POST) and (len(request.POST['email']) > 0):
            try:
                validate_email(request.POST['email'])
                email = request.POST['email']
            except:
                return {'status':'error', 'type':'email', 'message': 'Invalid email'}
        else:
            if request.user.is_authenticated():
                email = request.user.email
            else:
                email = ""
        if ("job-comments" in request.POST):
            if (len(request.POST['job-comments']) == 0):
                job_name = ''
            else:
                if (len(request.POST['job-comments']) < 200):
                    job_name = request.POST['job-comments']
                else:
                    return {'status':'error', 'type':'job-comments', 'message': 'Comments must be &lt;200 characters'}
        else:
            job_name = ""
        if ('email-subject' in request.POST):
            if (len(request.POST['email-subject']) == 0):
                email_subject = ''
            else:
                if (len(request.POST['email-subject']) < 100):
                    email_subject = request.POST['email-subject']
                else:
                    return {'status':'error', 'type':'email-subject', 'message': 'No more than 100 characters'}
        params = {}
        try:
            params['stage1-eval'] = float(request.POST['stage-one-eval'])
            params['stage1-mda-family-qcov'] = float(request.POST['stage-one-mda-qcov'])/100
            params['stage1-mda-family-hcov'] = float(request.POST['stage-one-mda-hcov'])/100
            params['stage1-pfam-family-hcov'] = float(request.POST['stage-one-pfam-hcov'])/100
            params['stage2-eval'] = float(request.POST['stage-two-eval'])
            params['stage2-mda-family-qcov'] = float(request.POST['stage-two-mda-qcov'])/100
            params['stage2-mda-family-hcov'] = float(request.POST['stage-two-mda-hcov'])/100
            params['stage2-pfam-family-hcov'] = float(request.POST['stage-two-pfam-hcov'])/100
            params['query-hmm-pid'] = float(request.POST['query-hmm-pid'])
            params['ec-checkboxes'] = request.POST['ec-checkboxes'].strip()
            params['ec-kerf-threshold'] = int(request.POST['ec-kerf-threshold'])
            params['ec-orthology-number'] = int(request.POST['ec-orthology-number'])
            params['cluster-similarity'] = int(request.POST['cluster-similarity'])
            params['stage2-maxseqs'] = int(request.POST['stage-two-maxseqs'])
            params['ortholog-coverage'] = int(request.POST['ortholog-coverage'])
            params['query-coverage'] = int(request.POST['query-coverage'])
            params['minimum-pid-for-orthology'] = int(request.POST['minimum-pid-for-orthology'])
            params['stage3-koa'] = int(request.POST['stage3-koa'])
            params['lambda'] = float(request.POST['consensus-uniprot-description-lambda'])
            params['amplitude'] = float(request.POST['consensus-uniprot-description-amplitude'])
            params['consensus-uniprot-threshold-high'] = float(request.POST['consensus-uniprot-threshold-high'])
            params['consensus-uniprot-threshold-medium'] = float(request.POST['consensus-uniprot-threshold-medium'])            
            ellis_island_hcov = float(request.POST['ellis-island-hcov'])/100
            ellis_island_qcov = float(request.POST['ellis-island-qcov'])/100
            ellis_island_maxseqs = int(request.POST['ellis-island-maxseqs'])
            ellis_island_maxfams = int(request.POST['ellis-island-maxfams'])
            ellis_island_taxon_pid = float(request.POST['ellis-island-taxon-pid'])/100
            
            ellis_island_string = str(request.POST['run-ellis-island'])
            if ellis_island_string.upper() == 'TRUE':
                run_ellis_island = True
            else:
                run_ellis_island = False
        except:
            return {'status': 'error', 'type': 'parameter', 'message': 'Malformed parameters.  Please check program parameters below and submit again.'}
        if params == HIGH_RECALL_PARAMS:
            parameters_chosen = 'High recall'
        elif params == HIGH_PRECISION_PARAMS:
            parameters_chosen = 'High precision'
        elif params == REMOTE_HOMOLOG_PARAMS:
            parameters_chosen = 'Remote homologs'
        elif params == FRAGMENT_PARAMS:
            parameters_chosen = 'Partial sequence'
        else:
            parameters_chosen = 'Custom'
        # This is a terrible fix for a documented problem with django.  We should upgrade django.
        # Until then, do this.
        # see https://groups.google.com/forum/?fromgroups=#!topic/django-developers/YbOHeu3s7vo
        try:
            job = Fatcat2Job.objects.get(id=1)
            job.fasta = 'A'
            job.save()
        except:
            pass
        try:
            oj = OhanaJobs.objects.get(id=1)
            oj.job_type='fatcat2'
            oj.save()
        except:
            pass
        job = Fatcat2Job.objects.create(
                fasta_sequence = fasta,
                status_id = 1,
                user_email = email,
                stage1_e_value = params['stage1-eval'],
                stage1_mda_family_query_coverage = params['stage1-mda-family-qcov'],
                stage1_mda_family_hmm_coverage = params['stage1-mda-family-hcov'],
                stage1_pfam_family_hmm_coverage = params['stage1-pfam-family-hcov'],
                stage2_e_value = params['stage2-eval'],
                stage2_mda_family_query_coverage = params['stage2-mda-family-qcov'],
                stage2_mda_family_hmm_coverage = params['stage2-mda-family-hcov'],
                stage2_pfam_family_hmm_coverage = params['stage2-pfam-family-hcov'],
                stage2_query_hmm_pid = params['query-hmm-pid'],
                stage2_orthology_methods = params['ec-checkboxes'],
                stage2_enclosing_clade_kerf_threshold = params['ec-kerf-threshold'],
                stage2_redundant_methods_required = params['ec-orthology-number'],
                stage2_max_sequences = params['stage2-maxseqs'],
                stage3_ortholog_cluster_similarity = params['cluster-similarity'],
                stage3_ortholog_coverage = params['ortholog-coverage'],
                stage3_query_coverage = params['query-coverage'],
                stage3_query_ortholog_pid = params['minimum-pid-for-orthology'],
                stage3_max_sequences = params['stage3-koa'],
                stage4_consensus_uniprot_weighting_lambda = params['lambda'],
                stage4_consensus_uniprot_weighting_amplitude = params['amplitude'],
                stage4_consensus_uniprot_threshold_for_high_confidence = params['consensus-uniprot-threshold-high'],
                stage4_consensus_uniprot_threshold_for_medium_confidence = params['consensus-uniprot-threshold-medium'],
                fatcat_job_name = job_name,
                email_subject = email_subject,
                save_job = False,
                parameterization = parameters_chosen,
                is_done = False,
                test = False,
                run_ellis_island = run_ellis_island,
                ellis_island_maxseqs = ellis_island_maxseqs,
                ellis_island_maxfams = ellis_island_maxfams,
                ellis_island_member_coverage = ellis_island_hcov,
                ellis_island_query_coverage = ellis_island_qcov,
                ellis_island_taxon_clustering_pid = ellis_island_taxon_pid
        )
        if seq_header:
            header = seq_header
        else:
            header = "FAT-CAT job %d" % job.id
        job.fasta_header = header
        job.save()
        # this isn't completely thought out...
        # log this job in the job table.
        ojob = OhanaJobs.objects.create(
            submitting_user = u,
            submitting_ip_address = ip,
            job_type = 'fatcat2'
        )
        # reference the ohana job entry in the fatcat job table
        job.ohana_job = ojob
        job.save()
        return self.read(None, job.id)

class Fatcat2Functions(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        return

# all the functions, no fancy stuff
class Fatcat2AllFunctions(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        return {'functions': job.functions}
        '''iTotalRecords = job.functions.filter(cluster_representative__classification='Orthologs').count()        

        #if sortingDirection == 'desc':
        #    sort_string = '-' + sort_string

        slice = job.clusters.filter(cluster_representative__classification='Orthologs').order_by('cluster_representative__member_to_query_pid')
        
        fatcat_clusters = list(slice)

        aaData = []

        for cluster in fatcat_clusters:
            aaData.append(cluster.get_candidate_ortholog_row_datatable_object())
        pprint.pprint(aaData)
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords}'''
        
class Fatcat2OtherSequenceMatches(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        sEcho = int(request.GET.get('sEcho'),0)
        iDisplayLength = int(request.GET.get('iDisplayLength', 10))
        startRecord = int(request.GET.get('iDisplayStart', 0))
        endRecord = startRecord + iDisplayLength
        numCols = int(request.GET.get('iSortingCols',0))
        sortedColumn = None
        # get the column we want to sort on - now this only sorts by one column
        # can add more sorting columns later.
        sortedColumn = int(request.GET.get('iSortCol_0', 0))
        sortingDirection = request.GET.get('sSortDir_0', 'asc')
        iTotalRecords = iTotalDisplayRecords = job.clusters.filter(cluster_representative__classification='Other Sequence Matches').count()        

        sort_string = None

        # which fields?
        if sortedColumn == 0:
            # sort by family type
            sort_string = 'cluster_representative__uniprot__in_swissprot_f'
        elif sortedColumn == 1:
            sort_string = 'cluster_representative__uniprot__uniprot_identifier'
        elif sortedColumn == 2:
            sort_string = 'cluster_representative__uniprot__description'
        elif sortedColumn == 3:
            sort_string = 'cluster_representative__uniprot__uniprotgene__name'
        elif sortedColumn == 4:
            sort_string = 'taxon__scientific_name'
        elif sortedColumn == 5:
            sort_string = 'cluster_representative__member_to_query_pid'
        elif sortedColumn == 6:
            sort_string = 'cluster_representative__query_coverage'
        elif sortedColumn == 7:
            sort_string = 'cluster_representative__member_coverage'
        else:
            return {}
           
        if sortingDirection == 'desc':
            sort_string = '-' + sort_string

        slice = job.clusters.filter(cluster_representative__classification='Other Sequence Matches').order_by(sort_string)[startRecord:endRecord]
        
        fatcat_clusters = list(slice)

        aaData = []

        for cluster in fatcat_clusters:
            aaData.append(cluster.get_candidate_ortholog_row_datatable_object())
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords, 'iTotalDisplayRecords': iTotalDisplayRecords, 'sEcho': sEcho}

class Fatcat2CandidateOrthologs(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        sEcho = int(request.GET.get('sEcho'),0)
        iDisplayLength = int(request.GET.get('iDisplayLength', 10))
        startRecord = int(request.GET.get('iDisplayStart', 0))
        endRecord = startRecord + iDisplayLength
        numCols = int(request.GET.get('iSortingCols',0))
        sortedColumn = None
        # get the column we want to sort on - now this only sorts by one column
        # can add more sorting columns later.
        sortedColumn = int(request.GET.get('iSortCol_0', 0))
        sortingDirection = request.GET.get('sSortDir_0', 'asc')
        iTotalRecords = iTotalDisplayRecords = job.clusters.filter(cluster_representative__classification='Orthologs').count()        

        sort_string = None

        # which fields?
        if sortedColumn == 0:
            # sort by family type
            sort_string = 'cluster_representative__uniprot__in_swissprot_f'
        elif sortedColumn == 1:
            sort_string = 'cluster_representative__uniprot__uniprot_identifier'
        elif sortedColumn == 2:
            sort_string = 'cluster_representative__uniprot__description'
        elif sortedColumn == 3:
            sort_string = 'cluster_representative__uniprot__uniprotgene__name'
        elif sortedColumn == 4:
            sort_string = 'taxon__scientific_name'
        elif sortedColumn == 5:
            sort_string = 'cluster_representative__member_to_query_pid'
        elif sortedColumn == 6:
            sort_string = 'cluster_representative__query_coverage'
        elif sortedColumn == 7:
            sort_string = 'cluster_representative__member_coverage'
        else:
            return {}
           
        if sortingDirection == 'desc':
            sort_string = '-' + sort_string

        slice = job.clusters.filter(cluster_representative__classification='Orthologs').order_by(sort_string)[startRecord:endRecord]
        
        fatcat_clusters = list(slice)

        aaData = []

        for cluster in fatcat_clusters:
            aaData.append(cluster.get_candidate_ortholog_row_datatable_object())
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords, 'iTotalDisplayRecords': iTotalDisplayRecords, 'sEcho': sEcho}

# returns all the orthologs, no pagination etc. 
class Fatcat2AllCandidateOrthologs(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        iTotalRecords = job.clusters.filter(cluster_representative__classification='Orthologs').count()        

        #if sortingDirection == 'desc':
        #    sort_string = '-' + sort_string

        slice = job.clusters.filter(cluster_representative__classification='Orthologs').order_by('cluster_representative__member_to_query_pid')
        
        fatcat_clusters = list(slice)

        aaData = []

        for cluster in fatcat_clusters:
            aaData.append(cluster.get_candidate_ortholog_row_datatable_object())
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords}
        
class Fatcat2EnclosingClades(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        sEcho = int(request.GET.get('sEcho'),0)
        iDisplayLength = 10 # we are only doing 10 now...
        startRecord = int(request.GET.get('iDisplayStart', 0))
        endRecord = startRecord + iDisplayLength
        numCols = int(request.GET.get('iSortingCols',0))
        sortedColumn = None
        # get the column we want to sort on - now this only sorts by one column
        # can add more sorting columns later.
        sortedColumn = int(request.GET.get('iSortCol_0', 0))
        sortingDirection = request.GET.get('sSortDir_0', 'asc')
        iTotalRecords = iTotalDisplayRecords = job.families.filter(passed_stage2=True).count()        

        sort_string = None

        # which fields?
        if sortedColumn == 0:
            # sort by the bpg accession
            sort_string = 'family__id'
        elif sortedColumn == 1:
            # sort by family type
            sort_string = 'family_type'
        elif sortedColumn == 2:
            sort_string = 'enclosing_clade_description'
        elif sortedColumn == 3:
            sort_string = 'enclosing_clade_taxonomy__scientific_name'
        elif sortedColumn == 4:
            sort_string = 'top_scoring_node_e_value'
        elif sortedColumn == 5:
            sort_string = 'query_hmm_consensus_pid'
        elif sortedColumn == 6:
            sort_string = 'match_data__stage2_query_coverage'
        elif sortedColumn == 7:
            sort_string = 'match_data__stage2_hmm_coverage'
        elif sortedColumn == 8:
            sort_string = 'match_data__stage2_region_aligned'
        else:
            return {}
           
        if sortingDirection == 'desc':
            sort_string = '-' + sort_string

        slice = job.families.filter(passed_stage2=True).order_by(sort_string)[startRecord:endRecord]
        
        fatcat_families = list(slice)

        aaData = []
        heatmap_data = []

        for family in fatcat_families:
            md = family.match_data.get()
            # append data for the heat map for these families
            heatmap_data.append({
                'id': family.family.id,
                'family_name': family.family_description,
                'node_name': family.top_scoring_node_description,
                'alignment_match_from': md.stage2_ali_from,
                'alignment_match_to': md.stage2_ali_to,
                'evalue': family.top_scoring_node_e_value,
                'node_id': family.top_scoring_node.id
            })
            # append data for the datatables record                
            aaData.append(family.get_enclosing_clade_row_datatable_object())
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords, 'iTotalDisplayRecords': iTotalDisplayRecords, 'sEcho': sEcho, 'heatmapData': heatmap_data}

class Fatcat2DistantClades(BaseHandler):
    alowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        sEcho = int(request.GET.get('sEcho'),0)
        iDisplayLength = 10 # we are only doing 10 now...
        startRecord = int(request.GET.get('iDisplayStart', 0))
        endRecord = startRecord + iDisplayLength
        numCols = int(request.GET.get('iSortingCols',0))
        sortedColumn = None
        # get the column we want to sort on - now this only sorts by one column
        # can add more sorting columns later.
        sortedColumn = int(request.GET.get('iSortCol_0', 0))
        sortingDirection = request.GET.get('sSortDir_0', 'asc')
        iTotalRecords = iTotalDisplayRecords = job.families.filter(passed_stage1=True, passed_stage2=False).count()        

        sort_string = None

        # which fields?
        if sortedColumn == 0:
            sort_string = 'reason_eliminated'
        elif sortedColumn == 1:
            # sort by the bpg accession
            sort_string = 'family__id'
        elif sortedColumn == 2:
            # sort by family type
            sort_string = 'family_type'
        elif sortedColumn == 3:
            sort_string = 'top_scoring_node_description'
        elif sortedColumn == 4:
            sort_string = 'top_scoring_node_taxonomy__scientific_name'
        elif sortedColumn == 5:
            sort_string = 'top_scoring_node_e_value'
        elif sortedColumn == 6:
            sort_string = 'query_hmm_consensus_pid'
        elif sortedColumn == 7:
            sort_string = 'match_data__stage2_query_coverage'
        elif sortedColumn == 8:
            sort_string = 'match_data__stage2_hmm_coverage'
        elif sortedColumn == 9:
            sort_string = 'match_data__stage2_region_aligned'
        else:
            return {}
           
        if sortingDirection == 'desc':
            sort_string = '-' + sort_string

        slice = job.families.filter(passed_stage1=True, passed_stage2=False).order_by(sort_string)[startRecord:endRecord]
        
        fatcat_families = list(slice)

        aaData = []
        heatmap_data = []

        for family in fatcat_families:
            md = family.match_data.get()
            # append data for the heat map for these families
            heatmap_data.append({
                'id': family.family.id,
                'family_name': family.family_description,
                'node_name': family.top_scoring_node_description,
                'alignment_match_from': md.stage2_ali_from,
                'alignment_match_to': md.stage2_ali_to,
                'evalue': family.top_scoring_node_e_value,
                'node_id': family.top_scoring_node.id
            })
            # append data for the datatables record                
            try:
                aaData.append(family.get_distant_clade_row_datatable_object())
            except:
                print family.id
                
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords, 'iTotalDisplayRecords': iTotalDisplayRecords, 'sEcho': sEcho, 'heatmapData': heatmap_data}

class Fatcat2FamilyMatches(BaseHandler):
    alowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}
        sEcho = int(request.GET.get('sEcho'),0)
        iDisplayLength = 10 # we are only doing 10 now...
        startRecord = int(request.GET.get('iDisplayStart', 0))
        endRecord = startRecord + iDisplayLength
        numCols = int(request.GET.get('iSortingCols',0))
        sortedColumn = None
        # get the column we want to sort on - now this only sorts by one column
        # can add more sorting columns later.
        sortedColumn = int(request.GET.get('iSortCol_0', 0))
        sortingDirection = request.GET.get('sSortDir_0', 'asc')
        iTotalRecords = iTotalDisplayRecords = job.families.all().count()        

        sort_string = None

        # which fields?
        if sortedColumn == 0:
            # sort by the bpg accession
            sort_string = 'family__id'
        elif sortedColumn == 1:
            # sort by family type
            sort_string = 'family_type'
        elif sortedColumn == 2:
            sort_string = 'family_description'
        elif sortedColumn == 3:
            sort_string = 'family_taxonomy__scientific_name'
        elif sortedColumn == 4:
            sort_string = 'family_e_value'
        elif sortedColumn == 5:
            sort_string = 'match_data__stage1_region_aligned'
        elif sortedColumn == 6:
            sort_string = 'match_data__stage1_hmm_coverage'
        elif sortedColumn == 7:
            sort_string = 'match_data__stage1_query_coverage'
        else:
            return {}
           
        if sortingDirection == 'desc':
            sort_string = '-' + sort_string

        slice = job.families.all().order_by(sort_string)[startRecord:endRecord]
        
        fatcat_families = list(slice)

        aaData = []
        heatmap_data = []

        for family in fatcat_families:
            md = family.match_data.get()
            # append data for the heat map for these families
            heatmap_data.append({
                'id': family.family.id,
                'family_type': family.family_type,
                'alignment_match_from': md.stage1_ali_from,
                'family_name': family.family_description,
                'alignment_match_to': md.stage1_ali_to,
                'evalue': family.family_e_value
            })
            # append data for the datatables record                
            aaData.append(family.get_family_row_datatable_object())
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords, 'iTotalDisplayRecords': iTotalDisplayRecords, 'sEcho': sEcho, 'heatmapData': heatmap_data}


class Fatcat2Cluster(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id, cluster_id):
        try:
            j = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}

        cluster_members = j.members.filter(cluster = int(cluster_id))
        if cluster_members[0].uniprot.taxon.common_name:
            common_name = cluster_members[0].uniprot.taxon.common_name
        else:
            common_name = ''

        # For now we are compiling the html in here, this is bad and makes me feel
        # really bad, but its so much harder not to do it this way.
        return_html = '<div class="cluster-detail"><h4>Members of cluster %d in <span' % (int(cluster_id))
        if common_name:
            return_html += ' class="tip" title="%s">%s</span></h4><table class="table"><thead><tr><th>SP</th><th>Uniprot ID</th><th>Description</th><th>Gene Name</th><th>%% ID</th><th>Q-Cov %%</th><th>H-Cov %%</th></tr></thead><tbody>' % (common_name, cluster_members[0].uniprot.taxon.scientific_name)
        else:
            return_html += '>%s</span></h4><table class="table"><thead><tr><th>SP</th><th>Uniprot ID</th><th>Description</th><th>Gene Name</th><th>%% ID</th><th>Q-Cov %%</th><th>H-Cov %%</th></tr></thead><tbody>' % (cluster_members[0].uniprot.taxon.scientific_name)

        for member in cluster_members:
            # get the gene name
            gene_name = UniProtGene.objects.filter(uniprot = member.uniprot)
            try:
                gene_name = gene_name[0].name
            except:
                gene_name = ''

            return_html += '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.1f</td><td>%.1f</td></tr>' % (member.get_swissprot_link(), member.get_uniprot_id_link(), member.uniprot.description, gene_name, member.get_member_pid_link(), member.query_coverage, member.member_coverage)
        return_html += '</tbody></div>'
        return {'table': return_html, 'cluster': cluster_id}

class Fatcat2GOAnnotation(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, fatcat_id, annotation_id):
        try:
            j = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}

        go_annotations_queryset = j.go.filter(annotation__ontology_term__accession='GO:'+annotation_id)
        sEcho = int(request.GET.get('sEcho'),0)
        iDisplayLength = int(request.GET.get('iDisplayLength', 10)) # we are only doing 10 now...
        startRecord = int(request.GET.get('iDisplayStart', 0))
        endRecord = startRecord + iDisplayLength
        numCols = int(request.GET.get('iSortingCols',0))
        sortedColumn = None
        # get the column we want to sort on - now this only sorts by one column
        # can add more sorting columns later.
        sortedColumn = int(request.GET.get('iSortCol_0', 0))
        sortingDirection = request.GET.get('sSortDir_0', 'asc')
        iTotalRecords = iTotalDisplayRecords = go_annotations_queryset.count()        
        sort_string = None

        # which fields?
        if sortedColumn == 0:
            # sort by the uniprot identifier
            sort_string = 'member__uniprot__uniprot_identifier'
        elif sortedColumn == 1:
            # sort by description
            sort_string = 'member__uniprot__description'
        elif sortedColumn == 2:
            sort_string = 'member__uniprot__taxon__scientific_name'
        elif sortedColumn == 3:
            sort_string = 'annotation__evidence__priority'
        else:
            return {}
           
        if sortingDirection == 'desc':
            sort_string = '-' + sort_string

        slice = go_annotations_queryset.order_by(sort_string)[startRecord:endRecord]
 
        these_go_annotations = list(slice)

        aaData = []

        for annotation in these_go_annotations:
            # append data for the datatables record                
            aaData.append(annotation.get_summary_row_datatable_object())
 
        return {'aaData': aaData, 'iTotalRecords':iTotalRecords, 'iTotalDisplayRecords': iTotalDisplayRecords, 'sEcho': sEcho}        

class Fatcat2Tree(BaseHandler):
    alowed_methods = ('GET')

    def read(self, request, fatcat_id):
        try:
            job = Fatcat2Job.objects.get(id=fatcat_id)
        except:
            return {}

        return {'tree': job.tree_string}        
