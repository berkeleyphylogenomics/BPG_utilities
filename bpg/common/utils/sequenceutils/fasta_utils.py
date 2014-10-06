#!/usr/bin/env python
'''
Basic functions for working with fasta files.
TODO: More sources for fasta files need to be added here.
'''
import re
import pfacts003

def get_uniprot_accession_from_header(header):
    '''Get uniprot accession from header regardless of source.'''
    if (header):
        try:
            inf = header.split()[0].lstrip('>')
        except:
            return None
        if '|' in inf:
            try:
                acc = inf.split('|')[1]
            except:
                return None
        elif ':' in inf:
            try:
                acc = inf.split(':')[1]
            except:
                return None
        elif '/' in inf:
            matches = re.search(r'(.*)/.*', inf)
            if matches:
                acc = matches.group(1)
            else:
                return None
        # this should be last, because it will match all uniprot identifiers...
        elif '_' in inf:
            try:
                acc = inf.split('_')[1]
            except:
                return None
        return acc
#        try:
#            j = pfacts003.phylofacts.models.UniProt.find_by_accession_or_identifier(acc)
#            return acc
#        except:
#            return None
    else:
        return None

def get_uniprot_from_header(header):
    '''Get uniprot object from header regardless of source.'''
    if (header):
        try:
            inf = header.split()[0].lstrip('>')
        except:
            return None
        if '|' in inf:
            try:
                acc = inf.split('|')[1]
            except:
                return None
        elif ':' in inf:
            try:
                acc = inf.split(':')[1]
            except:
                return None
        elif '/' in inf:
            matches = re.search(r'(.*)/.*', inf)
            if matches:
                acc = matches.group(1)
            else:
                return None
        # this should be last, because it will match all uniprot identifiers...
        elif '_' in inf:
            try:
                acc = inf.split('_')[1]
            except:
                return None
        try:
            return pfacts003.phylofacts.models.UniProt.find_by_accession_or_identifier(acc)
        except:
            return None
    else:
        return None

def get_uniprot_accession_from_fasta_file(fasta_file):
    '''Return a list of uniprot accessions in the fasta file.'''
    headers = [line for line in open(fasta_file).readlines()
               if line.startswith('>')]
    return [get_uniprot_accession_from_header(header)
            for header in headers]


def lcl_replace_fasta_file(fasta_file):
    '''Avoid issues with fastacmd not working with the created blastable db.
    Just replacing the entire start chars (between > and the colon or pipe)
    with >lcl|'''
    file_list = open(fasta_file).readlines()
    final_list = []
    for line in file_list:
        if line.startswith('>'):
            if ':' in line.split()[0]:
                line = re.sub('^.*?:', '>lcl|', line)
            elif '|' in line.split()[0]:
                line = re.sub('^.*?\|', '>lcl|', line)
        final_list.append(line)
    return ''.join(final_list)

