#!/usr/bin/env python
'''Get pfam related information for sequences, families and phogs.'''

import os
import sys
#There are some environment issues, so the first import might not work.
from bpg.common.utils import sql_patterns as sp
    
cursor = sp.connect_to_server()

def get_families_from_pfam_accession(pfam_accession):
    '''Given a pfam domain name get a list of families containing this domain.

    The output is a dictionary with two keys PFAM and GHG'''
    all_families = get_pfam_families_from_pfam_accession(pfam_accession)
    all_families.update(get_ghg_families_from_pfam_accession(pfam_accession))
    return all_families


def get_pfam_families_from_pfam_accession(pfam_accession):
    '''Given a pfam accession get all pfam families using this accession as seed.'''
    pfam_matches = sp.sql_results(cursor,
                                  sp.PFAM_ACCESSION_PFAM_FAMILY,
                                  (pfam_accession,))
    pfam_matches = [pfam_match[0] for pfam_match in pfam_matches]
    return_dict = {'PFAM' : list(set(pfam_matches))}
    return return_dict


def get_ghg_families_from_pfam_accession(pfam_accession):
    '''Given a pfam accession get all GHG families containing this domain.'''
    ghg_matches = sp.sql_results(cursor,
                                  sp.PFAM_ACCESSION_GHG_FAMILY,
                                  (pfam_accession,))
    ghg_matches = [ghg_match[0] for ghg_match in ghg_matches]
    return_dict = {'GHG' : list(set(ghg_matches))}
    return return_dict


def get_pfam_accession_from_pfam_name(pfam_name):
    '''Get the pfam accession for a pfam name.'''
    pfam_accession = sp.sql_results(cursor,
                                    sp.PFAM_ACCESSION_FROM_NAME,
                                    (pfam_name,24))
    try:
        return pfam_accession[0][0]
    except:
        return None
    
def get_families_from_pfam_name(pfam_name):
    '''Get all families from pfam name.'''
    pfam_accession = get_pfam_accession_from_pfam_name(pfam_name)
    if pfam_accession:
        return get_families_from_pfam_accession(pfam_accession)
    else:
        return {'PFAM' : [],
                'GHG' : []}

if __name__ == "__main__":
    print get_families_from_pfam_accession('PF12417')
    print get_pfam_accession_from_pfam_name('CARD')
    print get_families_from_pfam_name('CARD')
