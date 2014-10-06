#!/usr/bin/python
'''
Utilities for getting phog and sequence relations.
'''

import sys




sys.path.append('/home/awarrier/unfuddle/phylofacts/ohana/TRUNK/bpg/common/utils')
import sql_patterns as sp
    
cursor = sp.connect_to_server()
def get_families_from_sequence(uniprot_accession):
    '''Get families containing given protein.'''
    return sp.sql_results(cursor,
                          sp.FAMILY_SEQUENCE_SQL,
                          (uniprot_accession,))

def get_sequences_for_families(family_accession_list):
    '''Get sequences in families.'''
    family_accession_as_tuple = '(%s)' % (','.join(family_accession_list))
    sequence_gathering_script = sp.SEQUENCE_FAMILY_SQL % family_accession_as_tuple
    return sp.sql_results_without_params(cursor,
                                         sequence_gathering_script)
                                         

    
def get_phogs_from_family(family):
    '''Get PHOGs for family id'''
    pass

def get_best_GHG_book_for_sequence(uniprot_accession):
    pass

def get_best_pfam_book_for_sequence(uniprot_accession):
    pass

if __name__ == "__main__":
    print get_families_from_sequence('O14727')
