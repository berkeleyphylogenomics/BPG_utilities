#!/usr/bin/env python2.7
'''Create a config file for a reaction. If there is a list_of_pf_families file in the reaction directory, then create a config file that will be used to populate the biocyc reaction pages.

This file is run assuming that the reaction_page_fam_list_collection.py script has already been run for the given reaction.

The code below is mostly identical to the code used for creating PhyloFacts-Pfam pages.
'''

import os
import sys
sys.path.insert(0,'/home/awarrier/ohana_repository')
from pfacts003.phylofacts.models import Family, MetacycReaction
from bpg.common.utils import sql_patterns as sp
cursor = sp.connect_to_server()

class Biocyc_WEB_DICT_CREATOR:
    def __init__(self, biocyc_process):
        self.biocyc_process = biocyc_process
        print 'Biocyc process: %s' % biocyc_process
        self.input_file = '/clusterfs/ohana/bpg/biocyc_rxn_pf_families/%s/list_of_pf_families' % biocyc_process
        self.output_file = '/clusterfs/ohana/bpg/biocyc_rxn_pf_families/%s/web_dict' % biocyc_process
        if os.path.exists(self.output_file):
            print 'The config file already exists for %s' % biocyc_process
            sys.exit()
        self.bpg_accessions = [line.strip() for line in open(self.input_file).readlines()]
        self.output_config_dictionary = {}
        self.output_config_dictionary['biocyc_process'] = biocyc_process
        metacyc_processes = MetacycReaction.objects.filter(id=biocyc_process)
        # Get a name for this reaction
        metacyc_common_name = None
        for process in metacyc_processes:
            if metacyc_common_name: break
            if process.common_name:
                metacyc_common_name = process.common_name
            else:
                metacyc_common_name = process.systematic_name
        self.output_config_dictionary['biocyc_common_name'] = metacyc_common_name
        self.output_config_dictionary['phylofacts_families'] = str(self.bpg_accessions)
        self.output_config_dictionary['family_info'] = {}
        self.output_config_dictionary['family_info']['values'] = []
        self.output_config_dictionary['num_unique_taxa'] = self.get_unique_taxa_within_family_list()
        self.output_config_dictionary['num_unique_sequences'] = self.get_unique_sequences_within_family_list()
        self.all_sequences = []
        
    def get_uniprot_descriptions_for_family(self):
        value_dict = self.output_config_dictionary['family_info']['values']
        for accession in self.bpg_accessions:
            print accession
            temp_dict = {}
            bpg_identifier = accession[-7:]
            temp_dict['bpg_accession'] = accession
            bpg_identifier = int(bpg_identifier)
            f = Family.find_by_family_id(bpg_identifier)
            description_tuples = f.canonical_root_node().get_description(returnAll=True, score_pfam=True)
            temp_dict['taxonomic_distribution'] = str(f.get_taxon())
            temp_dict['num_members'] = f.get_num_informative_members()
            temp_dict['num_unique_taxa'] = f.canonical_tree.number_of_taxa()
            temp_dict['uniprot_descriptions'] = '; '.join(self.get_top_five_from_descriptions(description_tuples)
                                                          )
            self.get_family_sequences(f)
            value_dict.append(temp_dict)

    def get_top_five_from_descriptions(self, description_tuples):
        '''Get the top scoring five.'''
        ret_list = self.get_all_descriptions(description_tuples)
        return ret_list[:10]

    def get_unique_taxa_within_family_list(self):
        '''Get unique taxa within the Biocyc process PhyloFacts families'''
        bpg_accessions = [accession[-7:] for accession in self.bpg_accessions]
        bpg_accessions = '(%s)' % ', '.join(bpg_accessions)
        unique_list = sp.sql_results_without_params(cursor, """select DISTINCT(sequence_header_taxon_id) from family_sequence_taxa_2012_02_17 where family_id in %s""" % bpg_accessions)
        return len([result[0] for result in unique_list])

    def get_unique_sequences_within_family_list(self):
        '''Get unique taxa within the Biocyc process PhyloFacts families'''
        bpg_accessions = [accession[-7:] for accession in self.bpg_accessions]
        bpg_accessions = '(%s)' % ', '.join(bpg_accessions)
        unique_list = sp.sql_results_without_params(cursor, """select DISTINCT(uniprot_accession) from family_sequence_taxa_2012_02_17 where family_id in %s""" % bpg_accessions)
        return len([result[0] for result in unique_list])
    
        
    def get_all_descriptions(self, description_tuples):
        '''Get all descriptions'''
        sortedlist = sorted(description_tuples, key=lambda score: score[1], reverse=True)
        retlist =[]
        for desc, score, raw in sortedlist:
            r = desc.strip()
            if r[-1] == ";":
                retlist.append((str(r[0:-1]), score, raw))
            else:
                retlist.append((str(r), score, raw))
        retlist = [tup[0] for tup in retlist]
        return retlist

    def get_family_sequences(self, family):
        '''Get family sequences as uniprot objects.'''
        leaves = family.canonical_root_node().get_included_leaves()
        self.all_sequences += [leaf.sequence_header.uniprot for leaf in leaves if leaf.sequence_header.uniprot is not None]

    def score_sequences(self):
        '''Score sequences in the collection of sequences'''
        raw_votes = {}
        votes_for_description = {}
        votes_for_informative_description = {}
        for description in self.all_sequences:
            if not raw_votes.has_key(description):
                raw_votes[description] = 1
            else:
                raw_votes[description] += 1
            if (uncharacterized_re.search(description) or predicted_re.search(description) or 
                hypothetical_re.search(description)):
              if not votes_for_uninformative_description.has_key(description):
                  votes_for_uninformative_description[description] = 0
              votes_for_uninformative_description[description] += 1
            else:
                if description not in votes_for_description:
                    votes_for_description[description] = 0
                if description.in_swissprot_f:
                    votes_for_description[description] += 10
                else:
                    votes_for_description[description] += 1
        if len(votes_for_description) > 0:
            used_votes_for_description = votes_for_description
        else:
            used_votes_for_description = votes_for_uninformative_description
        winning_votes = 0
        winning_description = ''
        for key, value in used_votes_for_description.items():
            returnlist.append((key, value, raw_votes[key]))

    def main(self):
        self.get_uniprot_descriptions_for_family()
        f = open(self.output_file, 'w')
        f.write(str(self.output_config_dictionary))
        f.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'The Biocyc process identifier is needed.'
    pwdc = Biocyc_WEB_DICT_CREATOR(sys.argv[1])
    pwdc.main()


