#!/usr/bin/env python2.7
'''This script takes a metacyc reaction id as input andd writes down the list of bpg families that include sequences
which participate in this reaction. This can be considered as the first step towards updating the biocyc family page.
The second step would be to run the reaction_page_config_file_creation.py script (which takes the same input). In theory, both these scripts can be combined into one.

So, given an input "RXN0-6734", the below code does the following:
1. Check if this reaction has already been processed. In other words, check if a folder called "RXN0-6734" exists in /clusterfs/ohana/bpg/biocyc_rxn_pf_families. If it does, look for list_of_pf_families file within. If it doesnt exist, continue.
2. Find all the biocyc proteins and uniprot accessions for this reaction.
3. Find all families that contain these accessions.
'''
import os
import sys
import time
sys.path.append('/home/awarrier/unfuddle/phylofacts/ohana/TRUNK/bpg/common/utils')
import sql_patterns as sp
from familyutils import sequence_family as sf

class ReactionPageCreator:
    '''Main class'''
    def __init__(self, reaction_id):
        self.reaction_id = reaction_id
        self.cursor = sp.connect_to_server()
        self.start_time = time.time()
        self.output_dir = os.path.join('/clusterfs/ohana/bpg/biocyc_rxn_pf_families',
                                        self.reaction_id)
        self.output_file = os.path.join(self.output_dir, 'list_of_pf_families')

    def get_membership_data(self):
        '''Returns uniprot sequences that participate in the metacyc reaction.'''
        # The below sql query is fundamentally similar to the recursive call stored in the models.py
        # only reversed
        sql_query = '''
        WITH RECURSIVE membership(group_id, member_id)
        AS (
        SELECT
        metacyc.metacyc_membership.metacyc_group_id,
        metacyc.metacyc_membership.metacyc_member_id
        FROM
        metacyc.metacyc_membership, metacyc.metacyc_reaction
        WHERE
        metacyc_reaction.id = %s AND
        metacyc_membership.metacyc_group_id = metacyc_reaction.id
        UNION ALL
        SELECT
        mm.metacyc_group_id,
        mm.metacyc_member_id
        FROM
        metacyc.metacyc_membership mm, membership m
        WHERE
        m.member_id = mm.metacyc_group_id
        )
        SELECT group_id, member_id
        FROM membership
        '''
        print '%s: Initial sql query started ...' % self.reaction_id
        results = sp.sql_results(self.cursor, sql_query, (self.reaction_id,))
        print '%s: This took %s ...' % (self.reaction_id, time.time() - self.start_time)
        self.start_time = time.time()
        return results

    def filter_to_biocyc_proteins(self, membership_data):
        '''Filters the membership data and returns protein.'''
        object_set = set()
        for member, group in membership_data:
            object_set.add(member)
            object_set.add(group)
        object_list = list(object_set)
        output_list = []
        # Return only those objectst which don't have these keywords
        print '%s: This took %s ...' % (self.reaction_id, time.time() - self.start_time)
        self.start_time = time.time()
        return [obj for obj in object_list
                if not (('RXN' in obj) or
                ('PWY' in obj) or
                ('CPLX' in obj))]

    def get_uniprot_accessions_from_biocyc(self,biocyc_proteins):
        '''Given a biocyc protein list, get all the uniprot accessions that correspond to it.
        '''
        biocyc_proteins = ["'%s'" % prot for prot in biocyc_proteins]
        biocyc_proteins_sql = "(%s)" % ', '.join(biocyc_proteins)
        results = sp.sql_results_without_params(self.cursor,
            '''
            select uniprot_accession from
            metacyc.metacyc_protein
            where id in %s and
            uniprot_accession is not NULL''' % biocyc_proteins_sql)
        result_lists = [val[0] for val in results]
        output_list = []
        for i in result_lists:
            output_list += i
        print '%s: This took %s ...' % (self.reaction_id, time.time() - self.start_time)
        self.start_time = time.time()
        #Remove duplicates before returning
        return list(set(output_list))

                           
    def get_uniprot_accessions_for_rxn(self):
        '''Returns uniprot accessions for a given reaction id.'''
        db_membership = self.get_membership_data()
        biocyc_proteins = self.filter_to_biocyc_proteins(db_membership)
        uniprot_accessions = self.get_uniprot_accessions_from_biocyc(biocyc_proteins)
        return uniprot_accessions

    def get_families_from_sequence_list(self, uniprot_accessions_list):
        '''From the list of uniprot accessions, get the list of families containing these sequences.'''
        pf_families = sf.get_families_from_sequence_list(uniprot_accessions_list)
        print '%s: This took %s ...' % (self.reaction_id, time.time() - self.start_time)
        self.start_time = time.time()
        if not pf_families:
            return []
        return ['bpg%s' % str(fam[0]).zfill(7) for fam in pf_families]

    def write_down_output_file(self, list_vals):
        f = open(self.output_file, 'w')
        f.write('\n'.join(list_vals))
        f.close()
        
    def exit_if_output_file_exists(self):
        '''If the list of phylofacts families already exists, then exit, else build a new directory.'''
        if os.path.exists(self.output_file):
            print 'The reaction id (%s) has already been processed.' % self.reaction_id
            sys.exit()
            
    def main(self):
        '''Main function.'''
        self.exit_if_output_file_exists()
        print '%s: Start ...' % self.reaction_id
        uniprot_accessions_for_rxn = self.get_uniprot_accessions_for_rxn()
        bpg_families = self.get_families_from_sequence_list(uniprot_accessions_for_rxn)
        if len(bpg_families) > 0:
            try:
                os.mkdir(self.output_dir)
            except:
                pass
            self.write_down_output_file(bpg_families)
        return bpg_families

if __name__ == "__main__":
    try:
        rxn_id = sys.argv[1]
    except:
        print 'The reaction id needs to be provided as an argument.'
        sys.exit()
    rpc = ReactionPageCreator(rxn_id)
    rpc.main()
    
