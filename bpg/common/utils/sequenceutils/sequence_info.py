import glob
import os
import subprocess
import sys
import time
import urllib2
import re
from Bio.SeqUtils import CheckSum
from pfacts003.phylofacts.models import *
#There are some environment issues, so the first import might not work.
sys.path.append('/home/yshen/ohana_repository/bpg/common/utils')
import sql_patterns as sp
from sequenceutils import fasta_utils as fu 

    
cursor = sp.connect_to_server()

def check_if_accessions_in_db(sequence_ID, genome_file=None):
    '''Check if uniprot accession has been added to PF3.0. Else retrieve from uniprot.
    
    Returns a list containing either the uniprot accession, those for sequences
    matching the seguid identical to the input sequence, or an empty list if everything
    fails.
    '''
    uniprot_matches = sp.sql_results(cursor,
                                     "select id from uniprot where accession=%s;",
                                     (sequence_ID,))
    if len(uniprot_matches) == 0:
        print 'No database entries for uniprot accession %s' % sequence_ID
        identical_accessions = \
get_accessions_for_proteins_not_in_db(sequence_ID, genome_file)
        if identical_accessions:
            print 'Identical accessions are %s' % ','.join(identical_accessions)
        return identical_accessions
    else:
        return [sequence_ID]        

def get_accessions_for_proteins_not_in_db(sequence_ID, genome_file=None):
    '''Get accessions not found in the database'''
    if genome_file:
        file_name=os.path.split(genome_file)[1]
        lcl_file = '%s.lcl' % file_name
        o = open(lcl_file,'w')
        o.write(fu.lcl_replace_fasta_file(genome_file))
        o.close()
        try:
            os.system ('formatdb -i %s -n %s -o T' %(lcl_file,lcl_file))
            seq = subprocess.Popen(['fastacmd', '-d', lcl_file,
                                    '-s', sequence_id],
                                   stdout = subprocess.PIPE)
            fasta_file = '%s.fasta' % sequence_id
            o=open(fasta_file,'w')
            o.write(seq.stdout.read())
            o.close() 
        except:
            log(logger,"Unable to get %s from %s\n" % (sequence_ID, genome_file))
        response=open(fasta_file,'r')
        record = SeqIO.parse(response, 'fasta').next()
        response.close()
        os.remove(fasta_file)
        os.remove(lcl_file)
        os.remove('formatdb.log')
        for filename in glob.glob('%s.p??' % lcl_file) :
            os.remove(filename)
        
        seguid = CheckSum.seguid(record.seq)
        uniprot_accessions = sp.sql_results(cursor,
                                            "select accession from uniprot where seguid=%s;",
                                            (seguid,))
        uniprot_accessions = [str(acc[0]) for acc in uniprot_accessions]
        return uniprot_accessions
    else:
        print "No fasta file provided"            

def get_tree_and_tree_node_ids_from_accessions(uniprot_accessions):
    '''Given a list containing uniprot accessions, get the tree and tree_node ids.'''
    tree_tn_list = []
    for accession in uniprot_accessions:
        tree_tn_list += sp.sql_results(cursor,
                                sp.PHOG_ACCESSION_SQL,
                                (accession,))
    return tree_tn_list


def get_phogs_from_sequence(uniprot_accession, phog_id_return=True):
    '''Returns phog information (in the form of tree id and tree_node id from uniprot accession.'''
    uniprot_accessions = check_if_accessions_in_db(uniprot_accession)
    if len(uniprot_accessions) == 0:
        return []
    phog_tree_tn_list = get_tree_and_tree_node_ids_from_accessions(uniprot_accessions)
    if not phog_id_return:
        return phog_tree_tn_list
    else:
        if phog_tree_tn_list:
            return ['PHOG%s_%s' % (str(listvals[1]).zfill(7),
                                   str(listvals[2]).zfill(5))
                    for listvals in phog_tree_tn_list]

def get_GO_from_sequence_list(sequence_list):
    '''Returns a dictionary of GO terms given phog id

    Input should be a tuple string of sequences
    '''
    cursor.execute(sp.SEQUENCE_LIST_GO % sequence_list)
    go_results = cursor.fetchall()
    return go_results

def get_EC_from_sequence_list(sequence_list):
    '''Returns a list of EC numbers given a phog id.
    Input should be a tuple string of sequences
    '''
    cursor.execute(sp.SEQUENCE_LIST_EC % sequence_list)
    ec_results = cursor.fetchall()
    return ec_results

def get_KEGG_from_sequence_list(sequence_list):
    '''Returns a list of KEGG information given a phog id.
    Input can be either phog_id or a list [tree, left_id]
    '''
    ec_results = get_EC_from_sequence_list(sequence_list)
    ec_ids = [str(ec[0]) for ec in ec_results if not ec == []]
    if len(ec_ids) == 0:
        return []
    ec_ids = "(%s)" % ', '.join(ec_ids)
    print ec_ids
    cursor.execute(sp.EC_LIST_KEGG % ec_ids)
    kegg_results = cursor.fetchall()
    return kegg_results

def get_taxa_from_sequence_list(sequence_list):
    pass

def get_uniprot_descriptions_from_sequence_list(sequence_list):
    pass

if __name__ == "__main__":
    
    sequence_id = sys.argv[1]
    if len(sys.argv) <3:
      genome_file=None
    else:
      genome_file = sys.argv[2]
    print check_if_accessions_in_db(sequence_id, genome_file) 
