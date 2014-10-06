#!/usr/bin/python
'''Get the best pfam books for given Uniprot identifier from a genome.
Input argument is a uniprot identifier. Output is a tab delimited print of the following
Uniprot identifier, pfam domain name, pfam identifier, best phylofacts family, e_value

author:Ajith
The following steps need to be done before running this script:
1. Perform an hmmscan of the proteins from the selected genome into an accessible location.
e.g.
2. Create a blastable database from the reference proteome (so that any sequence can be obtained
using the identifier
The following steps are performed by the BestBookPFAM class (see the main function).
1.Creates a dictionary from hmmscan results in 1 above.
2.Gets pfam domains for given protein from this dictionary (names and accessions)
3.Gets all phylofacts 3.0 families for each domain
4.Gets the top scoring family (lowest e-value) when an hmmscan is run against hmms for each family.
5.Cleanup and remove all intermediate files except the final hmmscan output
'''
import os
import sys
import subprocess
import time
from bpg.common.parsers.hmm_parsers import parse_results_of_hmmsearch_or_hmmscan as hparser
from pfacts003.phylofacts.models import Family
# the location of bpg common utils might be different for you. Edit the below line accordingly.
sys.path.append('/home/awarrier/unfuddle/phylofacts/ohana/TRUNK/bpg/common/utils')
import sql_patterns as sp
import dir_of_family as dof
    
class BestBookPFAM:
    def __init__(self, uniprot_id, blastable_db):
        '''Initialize with the uniprot id.'''
        self.uniprot_id = uniprot_id
        self.blastable_db = blastable_db
        self.use_uniprot_pfams = True
        self.cursor = sp.connect_to_server()
        self.uniprot_sequence = self.get_uniprot_sequence(blastable_db)
        #1.
        try:
            os.mkdir('intermediates')
        except:
            pass
        #2.
        try:
            os.mkdir(os.path.join('intermediates', self.uniprot_id))
        except:
            pass
        #3.
        self.fasta_file = os.path.join('intermediates', self.uniprot_id, "%s.fa" % self.uniprot_id)
        f = open(self.fasta_file, 'w'); f.write(self.uniprot_sequence); f.close()

    def add_hmm_scan_and_blastdb(self, hmmscan_result):
        '''New function called only if hmmscan result is provided.
        1. Create an intermediate file directory if it does not exist
        2. Also create one for the temporary files for this protein. This will be removed
        in cleanup.
        3. Also write down the fasta sequence file for this uniprot id.'''
        self.hmmscan_result = hmmscan_result
                    
    def get_uniprot_sequence(self,blastable_db):
        '''Get sequence from uniprot id'''
        seq = subprocess.Popen(['fastacmd', '-d', blastable_db,
                                '-s', self.uniprot_id],
                               stdout = subprocess.PIPE)
        return seq.stdout.read()
        
    def generate_dict_of_hmmscan_result(self):
        '''If a dictionary of pfam scan results is not available, create it.

        To create, call hparser and get results (1). The dictionary has the unprot ids as
        keys and a list of pfam names as values. Write the string form of this dictionary (2).
        '''
        dict_stored_at = os.path.join('intermediates', 'hmmscan_dict')
        try:
            hmmscan_dict = open(dict_stored_at).read()
            hmmscan_dict = eval(hmmscan_dict)
            print 'Previous hmm dictionary found.'
        except:
            hmmscan_dict = {}
            #1.
            hmmscan_results = hparser.parse(self.hmmscan_result, 0.001, 1, 1)
            for query in hmmscan_results.hit_result_of_name_of_query:
                query_id = query.split(':')[-1]
                hmmscan_dict[query_id] = []
                for pfam_name in hmmscan_results.hit_result_of_name_of_query[query]:
                    hmmscan_dict[query_id].append(pfam_name)
            #2.
            f = open(dict_stored_at, 'w'); f.write(str(hmmscan_dict)); f.close()
            print 'New hmm dictionary created.'
        return hmmscan_dict

    def get_selected_pfam_domain_names(self, hmmscan_dict):
        '''Return the list of pfams associated with input uniprot id'''
        if not self.uniprot_id in hmmscan_dict:
            print 'Error:Absent in the hmm scan dictionary %s.' % self.uniprot_id
            sys.exit()
        return hmmscan_dict[self.uniprot_id]

    def get_selected_pfam_domain_names_from_uniprot(self):
        '''Return the list of pfams names associated with input uniprot id'''
        sql_search_string = """select pfam.name
        from uniprot_pfam, uniprot, pfam
        where uniprot.accession=%s and
        uniprot.id=uniprot_pfam.uniprot_id and
        uniprot_pfam.pfam_id = pfam.id;"""
        domain_matches = sp.sql_results(self.cursor, sql_search_string, (self.uniprot_id,))
        if domain_matches == []:
            return []
        else:
            return [domain[0] for domain in domain_matches]

    def convert_names_to_accessions(self, name_list):
        '''Returns the list of pfam accessions for the given list of pfam names.
        Check version 24 first(1). See the sql command.'''
        accessions = []
        for name_val in name_list:
            #1.
            accession_val = sp.sql_results(self.cursor, sp.PFAM_ACCESSION_FROM_NAME, (name_val, 24))
            if len(accession_val[0]) == 0:
                accession_val = sp.sql_results(self.cursor, sp.PFAM_ACCESSION_FROM_NAME, (name_val, 23))
            accessions.append(accession_val[0][0])
        return accessions

    def get_family_dict_for_pfam_accessions(self, accessions):
        '''Return a dictionary with input accessions as keys and corresponding PFAM families
        in a list as values. See sql query for details.
        If no families are found with any accession, stop here.'''
        family_dict = {}
        for accession in accessions:
            family_val = sp.sql_results(self.cursor, sp.PFAM_ACCESSION_PFAM_FAMILY, (accession,))
            family_list = [val[0] for val in family_val]
            if family_list == []:                
                continue
            family_dict[accession] = family_list
        if family_dict == {}:
            print 'None of the accessions had a family in %s.' % self.uniprot_id
            sys.exit()
        return family_dict

    def filter_family_dict_to_query_containing(self, family_dict):
        '''Filters the family_dict generated by the previous step.

        Filtering is done so that only bpg families containing the query sequence are
        kept.
        '''
        from bpg.common.utils.familyutils import sequence_family as sf
        sequence_families = sf.get_families_from_sequence(self.uniprot_id)
        sequence_families = [val[0] for val in sequence_families]
        new_dict = {}
        if sequence_families == []:
            return []
        for pfam_accession, family_list in family_dict.items():
            family_set = set(family_list)
            sequence_family_set = set(sequence_families)
            sequence_containing = family_set.intersection(sequence_family_set)
            sequence_containing_list = list(sequence_containing)
            new_dict[pfam_accession] = sequence_containing_list
        return new_dict
            
    
    def get_best_family_dict_with_hmmscans(self):
        '''Given a family dictionary, get the best family by performing hmmscans.
        See individual functions for more details.'''
        best_family_dict = {}
        for pfam_accession, family_list in self.family_dict.items():
            best_family_dict[pfam_accession] = self.get_highest_scoring_family_with_hmms(pfam_accession,
                                                                                         family_list)
        return best_family_dict

    def get_best_family_dict_with_scores(self):
        '''Given a family dictionary, get the best family by using scores stored in the database.
        '''
        best_family_dict = {}
        for pfam_accession, family_list in self.family_dict.items():
            best_family_dict[pfam_accession] = self.get_highest_scoring_family_with_scores(family_list)
        return best_family_dict


    def get_highest_scoring_family_with_hmms(self, pfam_accession, family_list):
        '''Gets highest scoring family by performing hmmscan on concatenated hmms.'''
        self.current_hmm_output = 'intermediates/%s/%s.hmm.out' % (self.uniprot_id,pfam_accession)
        concatenated_hmm = self.concatenate_hmms(pfam_accession, family_list)
        self.run_hmmpress(concatenated_hmm)
        (highest_scoring_family, highest_score) = self.score_hmm(concatenated_hmm, family_list)
        return highest_scoring_family, highest_score

    def get_highest_scoring_family_with_scores(self, family_list):
        '''Gets highest scoring family by querying the database for family scores.'''
        highest_scoring_family = ''
        highest_score = 0
        for family in family_list:
            family_obj = Family.objects.get(id=family)
            family_score = int(family_obj.get_score())
            print "Family: %s Score: %s" % (family, family_score)
            if (family_score > highest_score):
                highest_scoring_family = 'bpg%s' % str(family).zfill(7)
                highest_score = family_score
        return str(highest_scoring_family), str(highest_score)

    def score_hmm(self, concatenated_hmm, family_list):
        '''Score the uniprot sequence against the concatenated hmm library'''
        os.system('hmmscan -o %s --tblout %s %s %s' % (self.current_hmm_output + '.temp',
                                                       self.current_hmm_output,
                                                       concatenated_hmm,
                                                       self.fasta_file))
        return self.get_best_scoring_family_from_table(family_list)

    def get_best_scoring_family_from_table(self, family_list):
        '''Get the first result from the hmm scan output table, If none, then
        get the first value from the family_list and give it an evalue of 10'''
        lines = [line for line in open(self.current_hmm_output).readlines()
                 if not line.startswith('#')]
        if len(lines) == 0:
            return 'bpg%s' % str(family_list[0]).zfill(7), '10'
        line1 = lines[0].split()
        family_id = line1[0]
        e_value = line1[4]
        return family_id, e_value
    
    def run_hmmpress(self, concatenated_hmm):
        '''Run hmmpress on the concatenated hmm'''
        os.system('hmmpress -f %s > %s' % (concatenated_hmm,
                                           '%s.press.out' % concatenated_hmm))
        
    def concatenate_hmms(self, pfam_accession, family_list):
        '''Get the locations of hmms for all families in the list.
        Concatenate these into one hmm in the intermediates folder'''
        concatenated_hmm = 'intermediates/%s/%s.hmm' % (self.uniprot_id,pfam_accession)
        location_list = ['cat']
        location_list += [os.path.join(dof.get_dir_of_family_id(family_id),
                                      'bpg%s.hmm' % str(family_id).zfill(7))
                         for
                         family_id in family_list]
        location_list += ['>', concatenated_hmm]
        os_cmd = ' '.join(location_list)
        os.system(os_cmd)
        return concatenated_hmm
        
    def cleanup(self):
        '''Clean up all the files created and remove the cursor object.'''
        self.cursor.close()
        intermediate_location = 'intermediates/%s' % (self.uniprot_id)
        for file_name in os.listdir(intermediate_location):
            if not (file_name.endswith('.fa') or file_name.endswith('.hmm.out')):
                os.remove(os.path.join(intermediate_location, file_name))

    def main(self):
        '''Main function that glues together other functions'''
	print 'Starting process on uniprot id: %s' % self.uniprot_id
        if not self.use_uniprot_pfams:
            print 'Obtaining dictionary from hmmscan results.'
            hmmscan_dict = self.generate_dict_of_hmmscan_result()
            print 'Getting PFAM domain names for input protein.'
            pfam_domain_names = self.get_selected_pfam_domain_names(hmmscan_dict)
        else:
            print 'Getting PFAM domain names from uniprot data.'
            pfam_domain_names = self.get_selected_pfam_domain_names_from_uniprot()
        #Quit if no pfam domains are found for this protein.
        if pfam_domain_names == []:
            print 'No pfam domains found in %s.' % self.uniprot_id
            sys.exit()
	print pfam_domain_names
        print 'Getting accessions for PFAM domains.'
        pfam_domain_accessions = self.convert_names_to_accessions(pfam_domain_names)
	print pfam_domain_accessions
        print 'Getting phylofacts 3.0 PFAM families matching each domain.'
        self.family_dict = self.get_family_dict_for_pfam_accessions(pfam_domain_accessions)
	print self.family_dict
        print 'Filtering phylofacts 3.0 PFAM families to ones containing the query.'
        self.family_dict = self.filter_family_dict_to_query_containing(self.family_dict)
        print self.family_dict
        print 'Getting best phylofacts 3.0 PFAM family for each domain.'
        #Earlier choice was to get one based on hmm scoring
        #self.best_family_dict = self.get_best_family_dict_with_hmms()
        #New one uses scores stored in the database.
        self.best_family_dict = self.get_best_family_dict_with_scores()
	print self.best_family_dict
        print 'Removing intermediate files.'
        self.cleanup()
        #Print a tab-delimited string
        for ind, pfam_domain_name in enumerate(pfam_domain_names):
            output_string = []
            output_string.append('Best Book:')
            output_string.append(self.uniprot_id)
            output_string.append(pfam_domain_name)
            pf_accession  = pfam_domain_accessions[ind]
            output_string.append(pf_accession)
            #If one pfam domain does not have associated phylofacts families,
            # but others do, then take care of it.
            try:
                output_string.append(self.best_family_dict[pf_accession][0])
                output_string.append(self.best_family_dict[pf_accession][1])
                print '\t'.join(output_string)
            except:
                print ('No phylofacts families found with pfam accession %s for uniprot id %s.'
                       % (pf_accession, self.uniprot_id))

if __name__ == "__main__":
    start_time = time.time()
    try:
        uniprot_id, blastable_db = sys.argv[1].split()
    except:
        print 'Error:The uniprot id argument and blastable_db is needed'
        sys.exit()
    bbp = BestBookPFAM(uniprot_id, blastable_db)
    #To add pfam hmm scan and blastable database arguments, you'll need to change the script.
    # Basically add this and fix the remaining code.
    # bbp.add_hmm_scan_and_blastdb('location_of_pfam_hmmscan')
    bbp.main()
    print 'Elapsed time: %s' % (time.time() - start_time)
