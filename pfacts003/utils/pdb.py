import os
from pfacts003.utils.hmm import hmmbuild, hmmsearch, parse_domtblout

# Functions for dealing with pdb stuff.
PDB_FASTA = '/home/cyrus_afrasiabi/PDB/pdb_seqres.txt'

def get_pdb_hits(query, dir, options=['-E', '1E-3']):
    # build an hmm from the query
    (hmm, err) = hmmbuild(query, 'qhmm')
    
    hmm_file = open(os.path.join(dir, 'query.hmm'),'w')
    hmm_file.write(hmmbuild(query, 'qhmm')[0])
    hmm_file.close()

    # score the pdb database
    ((normal, domain), code) = hmmsearch(hmm_file.name, PDB_FASTA, options)

    # parse the domain table
    dom = parse_domtblout(domain)

    # figure out what to do with different structure types and also repeat regions along protein (e.g. 
    # APAF_HUMAN gives 962 hits because 

    return (normal, domain)
