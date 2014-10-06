#/!/usr/bin/env python2.7
"""
This file contains a base parser for the hmmer domtblout format.

All it does is goes through the domtblout file and creates a list, one element per line.
Each element of the list is a dictionary with keys shown below, and values populated from the file.
The user takes this list, and generates whatever data structure they want for their task.
"""

import sys
import string
import pprint

def parse_hmm_domtblout(file_path):
    file = open(file_path, 'r')
    file_lines = file.readlines()
    file.close()
    returnlist = []
    for line in file_lines:
        if line[0] != '#':
            fields = line.split()
            if len(fields) >= 23:
                currentdict = {} 
                currentdict['target_name'] = fields[0]
                currentdict['target_accession'] = fields[1]
                currentdict['target_length'] = int(fields[2])
                currentdict['query_name'] = fields[3]
                currentdict['query_accession'] = fields[4]
                currentdict['query_length'] = int(fields[5])
                currentdict['full_seq_evalue'] = float(fields[6])
                currentdict['full_seq_score'] = float(fields[7])
                currentdict['full_seq_bias'] = float(fields[8])
                currentdict['domain_number'] = int(fields[9])
                currentdict['number_of_domains_in_query'] = int(fields[10])
                currentdict['domain_c_evalue'] = float(fields[11])
                currentdict['domain_i_evalue'] = float(fields[12])
                currentdict['domain_score'] = float(fields[13])
                currentdict['domain_bias'] = float(fields[14])
                currentdict['hmm_coord_from'] = int(fields[15])
                currentdict['hmm_coord_to'] = int(fields[16])
                currentdict['alignment_from'] = int(fields[17])
                currentdict['alignment_to'] = int(fields[18])
                currentdict['envelope_from'] = int(fields[19])
                currentdict['envelope_to'] = int(fields[20])
                currentdict['acc'] = float(fields[21])
                currentdict['target_description'] = string.join(fields[22:])
                returnlist.append(currentdict)
    return returnlist

if __name__ == "__main__":
    try:
        hmm_scan_domtblout_file = sys.argv[1]
    except:
        print "HMM domtblout file is required."
        sys.exit()
    pprint.pprint(parse_hmm_domtblout(hmm_scan_domtblout_file)) 
    sys.exit()
