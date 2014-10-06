#!/usr/bin/python
'''
Utilities for getting phog and sequence relations.
'''
import os
import sys
import time
import urllib2
import re
from Bio.SeqUtils import CheckSum
from pfacts003.phylofacts.models import *
#There are some environment issues, so the first import might not work.
from bpg.common.utils import sql_patterns as sp
from bpg.common.utils.sequenceutils.sequence_info import *
    
cursor = sp.connect_to_server()


def get_sequences_from_phog(input_phog, phog_id=True):
    '''Returns a list of sequences within a given phog.

    Input can be either phog_id or a list [tree, left_id]
    '''
    if not phog_id:
        tree_id, left_id = input_phog
    else:
        input_phog = input_phog.replace('PHOG', '')
        tree_id, left_id = [int(val) for val in input_phog.split('_')]
    tree_node_id = sp.sql_results(cursor,
                                  "select distinct(tree_node.id) from tree_node where tree_node.tree_id=%s and tree_node.left_id=%s;",
                                  (tree_id, left_id))
    if tree_node_id:
        tree_node_id = tree_node_id[0][0]
    sequences = sp.sql_results(cursor,
                               sp.PHOG_DATA_SQL,
                               (tree_node_id, tree_id))
    sequences = [sequence[2] for sequence in sequences]
    return sequences

def get_sequence_list_as_tuple_string(input_phog, phog_id=True):
    '''Returns a string version of the tuple of sequence list.'''
    sequence_list = get_sequences_from_phog(input_phog, phog_id)
    sequence_list = ["'%s'" % str(seq) for seq in sequence_list]
    sequence_list = "(%s)" % ', '.join(sequence_list)
    return sequence_list

def filter_phog_list(phog_list, sequence):
    '''Return the item from phog_list which contains sequence.'''
    for phog in phog_list:
        sequence_list_in_phog = get_sequences_from_phog(phog)
        if sequence in sequence_list_in_phog:
            return phog
        
def get_GO_from_phog(input_phog, phog_id=True):
    '''Returns a dictionary of GO terms given phog id

    Input can be either phog_id or a list [tree, left_id]
    '''
    sequence_list = get_sequence_list_as_tuple_string(input_phog, phog_id)
    return get_GO_from_sequence_list(sequence_list)


def get_EC_from_phog(input_phog, phog_id=True):
    '''Returns a list of EC numbers given a phog id.
    Input can be either phog_id or a list [tree, left_id]
    '''
    sequence_list = get_sequence_list_as_tuple_string(input_phog, phog_id)
    return get_EC_from_sequence_list(sequence_list)

def get_KEGG_from_phog(input_phog, phog_id=True):
    '''Returns a list of KEGG information given a phog id.
    Input can be either phog_id or a list [tree, left_id]
    '''
    ec_results = get_EC_from_phog(input_phog, phog_id)
    ec_ids = [str(ec[0]) for ec in ec_results if not ec == []]
    if len(ec_ids) == 0:
        return []
    ec_ids = "(%s)" % ', '.join(ec_ids)
    cursor.execute(sp.EC_LIST_KEGG % ec_ids)
    kegg_results = cursor.fetchall()
    return kegg_results

if __name__ == "__main__":
    tree_tn_list = get_phogs_from_sequence('P77333')
    print tree_tn_list
    print get_sequences_from_phog([242758,291], phog_id=False)
    print get_sequences_from_phog('PHOG0242758_00291')
    print get_GO_from_phog('PHOG0242758_00291')
    print get_EC_from_phog('PHOG0242758_00291')
    print get_KEGG_from_phog('PHOG0242758_00291')
