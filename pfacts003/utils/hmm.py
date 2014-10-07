'''
    Lots of fun stuff for using hmm*

    It might all be obsolete if/when we move to the hmm daemon goodness

    Hopefully you find it less fugly than the other ones.

    Why is this not in BioPython???

    - Jonathan
'''

import subprocess, re
from Bio import AlignIO
from cStringIO import StringIO

'''Input: FASTA as a string

Returns: (tbl_out, dom_tbl_out)'''
# Why not kwargs? - we want to be able to pass in --incdomT and such, which is unpythonic, and therefore punishable by death
def hmmscan(fasta, hmm_file, options=['-E', '1E-3']):
    #TODO input validation on E
    #args = ["/clusterfs/ohana/home/ddineen/stuff/hmmer-3.0/src/hmmscan", 
    args = ["/clusterfs/ohana/software/bin/hmmscan", 
        "-o", "/dev/null", 
        "--notextw", 
        "--tblout", "/dev/stdout", 
        "--domtblout", "/dev/stderr",
        ] + options + [ 
        hmm_file,
        "/dev/stdin" ]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (process.communicate(fasta), process.returncode)
    
def hmmsearch(query_hmm, seq_file, options=['-E', '1E-3']):
    # Inputs - query_hmm is hmm file either as a string with the hmm
    # or a string to the file handle
    # seq_file is a path to the sequence database

    try:
        # is query_hmm a file?
        f = open(query_hmm,'r')
        args = ["/clusterfs/ohana/software/bin/hmm3search", 
            "-o", "/dev/null", 
            "--notextw", 
            "--tblout", "/dev/stdout", 
            "--domtblout", "/dev/stderr" 
            ] + options + [ 
            f.name, seq_file]
        process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return (process.communicate(), process.returncode)
    except:
        args = ["/clusterfs/ohana/software/bin/hmm3search", 
            "-o", "/dev/null", 
            "--notextw", 
            "--tblout", "/dev/stdout", 
            "--domtblout", "/dev/stderr" 
            ] + options + [ 
            '/dev/stdin', seq_file]
        process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return (process.communicate(query_hmm), process.returncode)

def hmmfetch(hmm, name, outfile):
    '''Extracted hmm is saved in the outfile.'''
    args = ["/clusterfs/ohana/software/bin/hmm3fetch", 
        "-o", outfile,
        hmm, 
        str(name)
        ]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.stdout.read()

def hmmalign(fasta, hmm, options=['--outformat','A2M']):
    '''Input arguments - an hmmfile and a fasta string
    Output - alignment is returned as a string.'''
    args = ["/clusterfs/ohana/software/bin/hmm3align"] + options + [hmm, "/dev/stdin"]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(fasta)[0]

def hmmemit(hmm):
    '''Input argument - hmm file
    Output - consensus sequence as a string.'''    
    args = ["/clusterfs/ohana/software/bin/hmm3emit",
            "-c",
            hmm
            ]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.stdout.read()

def prettyalign(a2mstring):
    '''Input argument - a2m string
    Output - pretty aligned a2m file as a string.'''
    output_alignment = []
    for line in a2mstring.split("\n"):
        if not line.startswith(">"):
            output_alignment.append(re.sub(r'[a-z.]','', line))
        else:
            output_alignment.append(line)
    return '\n'.join(output_alignment)

def hmmbuild(msa, name):
    args = ["/clusterfs/ohana/software/bin/hmmbuild", 
        "-o", "/dev/null",
        "--informat", "afa",
        "-n", str(name),
        "/dev/stdout", 
        "-",
        ]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(msa)

   
def jackhmmer(fasta, seqdb, E, N):
    args = ["/clusterfs/ohana/software/bin/jackhmmer",
            "-o", "/dev/null",
            "-E", str(E),
            "--incE", str(E),
            "-N", str(N),
            "--domE", str(E),
            "--tblout", "/dev/null",
            "--domtblout", "/dev/stderr",
            "-A", "/dev/stdout",
            "-", seqdb]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(fasta)

def hmmpress(file_name):
    args = ["/clusterfs/ohana/software/bin/hmmpress", file_name]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    return process.returncode

def gather_homologs_with_jackhmmer(fasta, **kwargs):
    """ Gather homologs with jackhmmer and return a msa of the results. """
    seq_db = kwargs.get("seqdb", "/clusterfs/ohana/external/UniProt/current/protein")
    E = kwargs.get("E", "0.000001")
    N = kwargs.get("N", "2")
    (msa, dom) = jackhmmer(fasta, seq_db, E, N)
    dom = None
    msa_file = StringIO(msa)
    msa = None
    return list(AlignIO.parse(msa_file, "stockholm"))

def pfamAscan(fasta, options=None):
    """ Scan the fasta against pfam A HMMs, and return a list of the parsed domtbl """
    ((normal, dom), ec) = hmmscan(fasta, "/clusterfs/ohana/external/pfam/Pfam26.0/Pfam-A.hmm", options)
   
    eligible_pfams = [pfam for pfam in parse_domtblout(dom) if pfam.this_domain_i_value < 1]

    def contains_overlapping_domains(pfam_list):
        ''' this function determines if there are any overlapping domains in this list '''
        for p in pfam_list:
            for q in pfam_list:
                if (p != q):
                    if (set(range(p.ali_from, p.ali_to+1)).intersection(set(range(q.ali_from, q.ali_to + 1)))):
                        return True
        return False
    
    def overlaps(pfam_list, p):
        ''' If this domain overlaps '''
        for pfam in pfam_list:
            if p != pfam:
                if (set(range(p.ali_from, p.ali_to+1)).intersection(set(range(pfam.ali_from, pfam.ali_to + 1)))):
                    return True
        return False

    def remove_overlapping_pfams(pfam_list, p):
        ''' returns a list with any overlapping domains in pfam_list to p removed '''
        retlist = []
        for pfam in pfam_list:
            if pfam != p:
                if not (set(range(p.ali_from, p.ali_to+1)).intersection(set(range(pfam.ali_from, pfam.ali_to + 1)))):
                    retlist.append(pfam)
            else:
                retlist.append(pfam)
        return retlist

    ''' pick the highest scoring domain by bit score '''

    while contains_overlapping_domains(eligible_pfams):
        eligible_pfams = sorted(eligible_pfams, key = lambda eval: eval.this_domain_i_value)
        for p in eligible_pfams:
            if overlaps(eligible_pfams, p):
                eligible_pfams = remove_overlapping_pfams(eligible_pfams, p)
                break

    return eligible_pfams 

def top_k_families(fasta, k, type="ghg", eval='1e-4', pfam_overlap=0.7, mda_overlap=0.7):
    if type == "ghg":
        ((normal, dom), ec) = hmmscan(fasta, "/clusterfs/ohana/bpg/pfacts/phylofacts3_GHG_090811.hmms", [ '-E', str(eval) ] )
    else:
        ((normal, dom), ec) = hmmscan(fasta, "/clusterfs/ohana/bpg/pfacts/phylofacts3_domain_090811.hmms", [ '-E', str(eval) ])

#    normal_results = parse_tblout(normal)
#    normal_results.sort(reverse=True)
#    dom_results = parse_domtblout(dom)
#    dom_results.sort(reverse=True)
#
#    our_results = normal_results[0:k]
#    for result in our_results:
#        result.doms = []
#
#    for dr in dom_results:
#        result = next((result for result in our_results if dr.target_name == result.target_name), None)
#        if result:
#            result.doms.append(dr)
#
#    return our_results

    dom_results = parse_domtblout(dom)
    dom_results.sort(reverse=True)

    seen = set()
    our_results = []
    for result in dom_results:
        if ( result.ali_to - result.ali_from ) < 50:
            continue
        if result.target_name not in seen:
            our_results.append(result)
            seen.add(result.target_name)

    return our_results[0:k]

''' Finds the top family for each domain, basically

Doing this generically is suprisingly complicated.  Its actually a quite interesting problem,
but one for another day.'''
def top_family_per_region(fasta, type="pfam"):
    if type == "ghg":
        ((normal, dom), ec) = hmmscan(fasta, "/clusterfs/ohana/bpg/pfacts/phylofacts3_GHG_090811.hmms")
    else:
        ((normal, dom), ec) = hmmscan(fasta, "/clusterfs/ohana/bpg/pfacts/phylofacts3_domain_090811.hmms", [ '--domE', '1e-4', '-E', '1e-3' ])

    dom_results = parse_domtblout(dom)
    dom_results.sort(reverse=True)

    # hard coded values are bad.
    #the correct way is to patch hmmscan
    # I am not doing that today.
    dom_results = [ d for d in dom_results if d.this_domain_i_value < 1e-5 ]

    best_nonoverlapping = []
    for dom in dom_results:
        if any(best_nonoverlapping_dom.ali_from <= dom.ali_to and best_nonoverlapping_dom.ali_to >= dom.ali_from for best_nonoverlapping_dom in best_nonoverlapping):
            continue
        best_nonoverlapping.append(dom)

    return best_nonoverlapping

def parse_tblout(tblout):
    output = []
    for line in tblout.split("\n"):
        if line and (line[0] != "#" and line.rstrip() != ''):
            output.append(NormalHmmerRow(line))

    return output
   
class HmmerRow:
    def __lt__(self, other):
        return self.sort_field < other.sort_field
    def __gt__(self, other):
        return self.sort_field > other.sort_field
    def __eq__(self, other):
        return self.sort_field == other.sort_field
    def __eq__(self, other):
        return self.sort_field != other.sort_field
    def __le__(self, other):
        return self.sort_field <= other.sort_field
    def __ge__(self, other):
        return self.sort_field >= other.sort_field
     
class NormalHmmerRow(HmmerRow):
    def __init__(self, row):
        (self.target_name,
        self.target_accession,
        self.query_name,
        self.query_accession,
        self.full_e_value,
        self.full_score,
        self.full_bias,
        self.best_domain_e_value,
        self.best_domain_score,
        self.best_domain_bias,
        self.exp,
        self.reg,
        self.clu,
        self.ovm,
        self.env,
        self.dom,
        self.rep,
        self.inc,
        self.description) = row.split()

        # Floatify
        for field in ['full_e_value', 'full_score', 'full_bias', 
            'best_domain_e_value', 'best_domain_score', 'best_domain_bias',
            'exp', 'reg', 'clu', 'ovm', 'env', 'dom', 'rep', 'inc']:
            setattr(self, field, float(getattr(self, field))) 

    @property
    def sort_field(self):
        return self.full_score


def parse_domtblout(domtblout):
    output = []
    for line in domtblout.split("\n"):
        if line and (line[0] != "#" and line.rstrip() != ''):
            output.append(DomHmmerRow(line))

    return output
        

class DomHmmerRow(HmmerRow):
    def __init__(self, row):
        splitrow = row.split()
        (self.target_name,
        self.target_accession,
        self.target_length,

        self.query_name,
        self.query_accession,
        self.query_length,

        self.full_e_value,
        self.full_score,
        self.full_bias,

        self.this_domain_number,
        self.this_domain_of,
        self.this_domain_c_value,
        self.this_domain_i_value,
        self.this_domain_score,
        self.this_domain_bias,

        self.hmm_from,
        self.hmm_to,

        self.ali_from,
        self.ali_to,

        
        self.env_from,
        self.env_to,

        self.acc) = splitrow[0:22]

        self.description = ' '.join(splitrow[22:])

        # Floatify
        for field in ['full_e_value', 'full_score', 'full_bias', 
            'this_domain_c_value', 'this_domain_i_value', 'this_domain_score', 'this_domain_bias', 'acc' ]:
            setattr(self, field, float(getattr(self, field))) 

        # Intify
        for field in ['hmm_from', 'hmm_to', 'ali_from', 'ali_to',
            'env_from', 'env_to']:
            setattr(self, field, int(getattr(self, field))) 

    @property
    def sort_field(self):
        return self.this_domain_score
