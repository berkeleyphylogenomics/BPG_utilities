#!/usr/bin/python
'''Parse hmm scan results and return all results in a dictionary.

TODO: Currently only some fields are returned. This will have to be fixed.
author: Ajith
Date:16 Nov 2011
'''
import os
import sys
import string

import parse_results_of_hmmsearch_or_hmmscan as hparser

class HmmToDict:
    def __init__(self, hmmscan_result_file, threshold = 0.001):
        '''Initialize hmm_dict and pfam fields to be populated.'''
        self.hmmscan_result_file = hmmscan_result_file
        self.threshold = threshold
        self.hmm_dict = {}
        self.hit_characteristics = ['best_bias', 'best_evalue', 'best_score',
                                    'description', 'exp', 'full_bias', 'full_evalue',
                                    'full_score', 'inclusion', 'n', 'name'
                                    ]
        #List with data for individual matches. Does not include unaligned_hit
        #which will be included later
        self.hit_match_info = ['seq_from', 'seq_to', 'hmm_from', 'hmm_to',
                               'i_evalue']

    def get_hmm_scan_results(self):
        '''Gets hmm scan results using hparser'''
        return hparser.parse(self.hmmscan_result_file,
                             self.threshold, 1, 1)

    def get_query_hits_from_hmm_results(self):
        '''Get query information from the scan results'''
        return self.get_hmm_scan_results().hit_result_of_name_of_query

    def populate_match_values(self, query, hit, matchDict):
        '''Populate self.hmm_dict with individual match related information.'''
        uppercase_translation = string.maketrans(string.lowercase, string.uppercase)
        dotdash = '.-'
        for match_number in matchDict:
            self.hmm_dict[query][hit]['matches'][match_number] = {}
            match_result = matchDict[match_number]
            unaligned_hit = match_result.aligned_hit.translate(
                uppercase_translation, dotdash)
            for match_info in self.hit_match_info:
                value = eval('match_result.%s' % match_info)
                self.hmm_dict[query][hit]['matches'][match_number][match_info] = value
                self.hmm_dict[query][hit]['matches'][match_number]['unaligned_hit'] = unaligned_hit
    
    def populate_hmm_dict_with_hit_info(self, query, hit, hit_items):
        '''Populate self.hmm_dict with hit related information for each query'''
        self.hmm_dict[query][hit] = {}
        for characteristic in self.hit_characteristics:
            value = eval('hit_items.%s' % characteristic)
            self.hmm_dict[query][hit][characteristic] = value
        self.hmm_dict[query][hit]['matches'] = {}
        self.populate_match_values(query, hit, hit_items.matches)

    def main(self):
        '''Main function.'''
        query_hits = self.get_query_hits_from_hmm_results()
        for query in query_hits:
            self.hmm_dict[query] = {}
            for hit, hit_items in query_hits[query].items():
                self.populate_hmm_dict_with_hit_info(query, hit, hit_items)
        return self.hmm_dict
                
        
def parse_hmm_to_dict(hmmscan_result_file, threshold = 0.001):
    '''Parse the hmm file and return a dictionary with its values'''
    htd = HmmToDict(hmmscan_result_file, threshold)
    return htd.main()

if __name__ == "__main__":
    try:
        hmm_scan_result_file = sys.argv[1]
    except:
        print 'HMM Scan file is required.'
        sys.exit()

    print parse_hmm_to_dict(hmm_scan_result_file)
