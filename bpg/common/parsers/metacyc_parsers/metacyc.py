#!/usr/bin/env python2.7
'''
This script deals with creating inputs for database creation for ecoli proteins from ecocyc.
'''
import os
import sys
import subprocess
import re
import StringIO
import time
import sqlite3
import psycopg2
import psycopg2.extras
import pprint

from Bio import SeqIO
from Bio.SeqUtils import CheckSum

sys.path.insert(0, '/home/awarrier/ohana_repository')
from bpg.common.utils import sql_patterns as sp

def split_by_tab(string):
    return string.split('\t')

cursor = sp.connect_to_server()
        
class MetacycParser:
    '''Main class for parsing metacyc files in a directory.
    '''
    def __init__(self, dat_dir):
        '''Creates empty lists, makes working directories, sets run name'''
        self.protein_file = None
        self.protein_list = []
        self.gene_list = []
        self.reaction_list = []
        self.pathway_list = []
        self.no_pathway_protein_list = []
        self.membership_list = []

        self.dat_dir = dat_dir
        self.base_dir = 'data'
        self.make_dir(self.base_dir)
        path_vars = self.dat_dir.split(os.sep)
        self.run_name = path_vars[-3]
        self.run_dir = os.path.join(self.base_dir, self.run_name)
        print 'Run name: ',self.run_name
        print 'Data directory: ', self.dat_dir

    def make_dir(self, dir_name):
        try:
            os.mkdir(dir_name)
        except:
            pass
        
    def create_and_conn_to_sqlite_db(self, db_name):
        self.db = db_name
        self.remove_db_if_present()
        self.conn = sqlite3.connect(self.db)
        self.conn.text_factory = str
        self.cursor = self.conn.cursor()
        
    def parse_metacyc_file(self, filename):
        '''Parsing a metacyc dat file and returning a dictionary with values.
        This function is from the milo lab at
        http://code.google.com/p/milo-lab/source/browse/trunk/src/pygibbs/metacyc.py'''
        metacyc_file = open(filename, 'r')
        curr_field = ""
        field_map = {}
        line = '#'; 
        line_counter = 0
        entry2fields_map = {}
        while (line):
            #line = str(metacyc_file.readline().rstrip())
            line = metacyc_file.readline().rstrip()
            line_counter += 1
            if (line.startswith(('#', '/')) and not line.startswith('//') ):
                continue

            field = line.split(' - ')[0]
            value = "".join(line.split(' - ')[1:])

            if (field == "//"):
                entry = field_map["UNIQUE-ID"]
                entry2fields_map[entry] = field_map
                field_map = {}
            else:
                if (field != ""):
                    curr_field = field
                if (curr_field in field_map):
                    field_map[curr_field] = field_map[curr_field] + "\t" + value
                else:
                    field_map[curr_field] = value
        metacyc_file.close()
        return entry2fields_map

    def get_feature_list(self, dat_dict):
        feature_list = set()
        for features in dat_dict.values():
            for feature in features.keys():
                feature_list.add(feature)
        return list(feature_list)

    def convert_dict_into_table(self, dat_dict, table_name):
        '''Converts the dictionary created by the parse function into a database
        table.
        '''
        feature_list = self.get_feature_list(dat_dict)
        # Removes uppercase and spaces in features
        feature_list_joined = [
            '_'.join(feature.lower().split()) for feature in feature_list]
        #Replaces non word characters in feature
        feature_list_joined = [re.sub(r'[^\w]','_', feature) for feature in
                               feature_list_joined]
        #Create table
        self.cursor.execute("create table %s (%s)" %
                            (table_name,
                             ','.join(feature_list_joined)))
        self.conn.commit()
        for key, values in dat_dict.items():
            temp_list = []
            for feature in feature_list:
                feature_value = values.get(feature, '')
                feature_value = feature_value.replace('"', '')
                temp_list.append('"%s"' % feature_value)
            self.cursor.execute(""" insert into %s values
            (%s)""" % (table_name,
                       ','.join(temp_list)))
            self.conn.commit()
        
    def remove_db_if_present(self):
        try:
            os.remove(self.db)
        except:
            pass

    def get_column_name(self, filename):
        '''Replaces dashes from the filename with underscores
        and returns this value to be used as the column name in the
        sqlite database.'''
        filename = os.path.basename(filename)
        column_name = os.path.splitext(filename)[0]
        column_name = re.sub(r'[^\w]','_', column_name)
        return column_name
        
    def close_connection(self):
        self.cursor.close()
        self.conn.close()

    def get_value_set(self, table, column):
        '''Returns sets of values and duplicates found in column in the query table.'''
        self.cursor.execute("select %s from %s" % (column, table))
        values = self.cursor.fetchall()
        value_set = set()
        duplicates = set()
        for value in values:
            value = value[0].split('\t')
            for individual_value in value:
                if individual_value.strip() == '':
                    continue
                if individual_value in value_set:
                    duplicates.add(individual_value)
                else:
                    value_set.add(individual_value)
        return value_set, duplicates


    def determine_all_relationships(self):
        '''For every column in every table in the database,
        figure out the relationship and print it.'''
        tables = self.get_tables()
        for table in tables:
            output_dict = {}
            for column in self.get_columns(table):
                relationship = self.db_test1(table, column)
                if not relationship in output_dict:
                    output_dict[relationship] = []
                output_dict[relationship].append(column)
            print table
            for relationship in output_dict:
                print '\t%s' % relationship
                for column in output_dict[relationship]:
                    print '\t\t%s' % column

    def get_combined_set_from_tables(self, table_column_dict):
        set1 = set()
        for table, column in table_column_dict.items():
            value_set, duplicates = self.get_value_set(table, column)
            set1 = set1 | value_set
        return set1

    def get_tables(self):
        '''Gets the names of tables in the sqlite database.'''
        self.cursor.execute("""SELECT name FROM sqlite_master
        WHERE type='table'
        ORDER BY name""")
        return [tup[0] for tup in self.cursor.fetchall()]

    def get_columns(self, table):
        self.cursor.execute("select * from %s" % table)
        return [tup[0] for tup in self.cursor.description]

    def check_if_required_files_present(self, dir_name):
        '''Reads the directory, makes sure the list of dat files are present
        in it and returns a list of path names.'''
        # Test if a file containing AA fasta sequences are present, else quit.
        # No point in following through if no AA sequences are available.
        # If there is a file that ends with .aa, then choose that as the source
        # fasta file, else choose the file called protseq.fsa
        return_with_error = True
        for filename in os.listdir(dir_name):
            if filename.endswith('.aa'):
                self.protein_file = os.path.join(dir_name, filename) 
                return_with_error = False
            elif not self.protein_file:
                if filename == 'protseq.fsa.pin':
                    self.protein_file = os.path.join(dir_name, 'protseq.fsa')
                    return_with_error = False
        if return_with_error:
            print 'Protein file absent ... Quitting ...'
            return None
        print 'Protein fasta file found ...'
        print 'Protein file: %s' % self.protein_file
        #Now check for the other dat files
        dat_files = ['genes.dat', 'pathways.dat',
                     'proteins.dat', 'reactions.dat', 
                     'enzrxns.dat']
        sqlite_dat_files = [os.path.join(dir_name, dat_file) for
                            dat_file in dat_files]
        for dat_file in sqlite_dat_files:
            if not os.path.exists(dat_file):
                print '%s absent ... Quitting ...' % dat_file
                sys.exit()
        print 'All dat files are found ...'
        return sqlite_dat_files
    
    def _test_db_1(self, table, column, verbose=False):
        '''Checks whether within a table, a particular column shows
        one-to-one, one-to-many or many-to-many relationship.'''
        self.cursor.execute("select %s from %s" % (column, table))
        values = self.cursor.fetchall()
        # multiple values separated by tabs
        multiple_val_rows = [value for value in values if
                             len(value[0].split('\t')) > 1]
        if verbose:
            print multiple_val_rows
        print (
            'The %s table has %s rows with multiple tab separated values in the %s column.' %
               (table,
                len(multiple_val_rows),
                column))
        value_set, duplicates = self.get_value_set(table, column)
        print 'There are %s duplicates of these entries in the %s column.' % (
            len(duplicates), column)
        if verbose:
            print 'Duplicates: ', duplicates
        # Check relationship
        # The logic here is that
        # if a column has no tab delimited values in any row, and there are no duplicates
        #         then this is a one to one relationship
        # if a column has tab delimited values in any row, but values are not found
        #         in any other row, then this is a one to many relationship.
        # If there are duplicates in the column, but no tab separated
        #         values in any row, then this is a many to one relationship.
        # if there are duplicates in the column and tab delimited values,
        #         then this is a many to many relationship
        tab_dl_num = len(multiple_val_rows)
        duplicate_num = len(duplicates)
        if tab_dl_num == 0 and duplicate_num == 0:
            return 'one-to-one'
        if tab_dl_num != 0 and duplicate_num == 0:
            return 'one-to-many'
        if tab_dl_num == 0 and duplicate_num != 0:
            return 'many-to-one'
        if tab_dl_num != 0 and duplicate_num != 0:            
            return 'many-to-many'

    def _test_db_2(self, table1, column1, table2, column2):
        '''Checks how many values in column 1 in table1 are found in
        column2 in table2.'''
        column1_values, column1_duplicates = self.get_value_set(table1, column1)
        column2_values, column2_duplicates = self.get_value_set(table2, column2)
        print 'Total number of values in %s(%s): %s' % (column1, table1, len(column1_values))
        print 'Total number of duplicates in %s(%s): %s' % (column1, table1, len(column1_duplicates))
        print 'Total number of values in %s(%s): %s' % (column2, table2, len(column2_values))
        print 'Total number of duplicates in %s(%s): %s' % (column2, table2, len(column2_duplicates))
        print 'Values in %s(%s) not found in %s(%s): %s' % (column1, table1,
                                                            column2, table2,
                                                            len(column1_values - column2_values))
        print 'Values in %s(%s) not found in %s(%s): %s' % (column2, table2,
                                                            column1, table1,
                                                            len(column2_values - column1_values))
        print 'Values shared by both %s(%s) and %s(%s) (intersection): %s' % (column2, table2,
                                                            column1, table1,
                                                            len(column2_values & column1_values))

    def _test_db_3(self, table_column_dict1, table_column_dict2):
        '''Check if values in table_column_dict1 and table_column_dict2
        overlap.'''
        set1 = self.get_combined_set_from_tables(table_column_dict1)
        f = open('gene_ids_from_dat', 'w'); f.write('\n'.join(list(set1))); f.close()
        set2 = self.get_combined_set_from_tables(table_column_dict2)
        print 'Total number of values in set1: %s' % len(set1)
        print 'Total number of values in set2: %s' % len(set2)
        print 'Values in set1 not found in set2: %s' % (len(set1 - set2))
        print list(set1-set2)[:10]
        print 'Values in set2 not found in set1: %s' % (len(set2 - set1))        


    def create_protein_file(self):
        '''Create the protein tab delimited file by reading the sqlite database and
        the uniprot accession mapping file.'''
        self.cursor.execute('select unique_id, common_name from proteins')
        output = []
        for row in self.cursor:
            temp_output = []
            #Add protein unique_id
            temp_output.append(row[0])
            self.protein_list.append(row[0])
            #common name and abbrev name
            c_names = row[1].split('\t')
            if c_names[0].strip() == '':
                temp_output.append('\N')
            else:
                temp_output.append(c_names[0])
            #identifier that states the type of node
            temp_output.append('protein')
            # db name
            temp_output.append(self.run_name)
            #uniprot accession
            uniprot_accessions = self.uniprot_mapping.get(row[0], '')
            if uniprot_accessions == '':
                uniprot_accessions = '\N'
            else:
                uniprot_accessions = ['"%s"' % uniprot_acc
                                      for uniprot_acc in uniprot_accessions]
                uniprot_accessions = '{%s}' % ','.join(uniprot_accessions)
            temp_output.append(uniprot_accessions)
            output.append('\t'.join(temp_output))
        self.write_and_close(os.path.join(self.run_dir,
                                          'metacyc_protein.tab_delimited'),
                             output)
            
    def create_gene_file(self):
        '''Creates the gene file and adds associations in the membership_file.'''
        try:
            self.cursor.execute('select unique_id, common_name, product from genes')
        except:
            self.cursor.execute('select unique_id, product from genes')
        values = self.cursor.fetchall()
        gene_set = set()
        for value in values:
            gene_name = value[0]
            # for cases with no common name
            if len(value) == 3:
                common_name = value[1].split('\t')[0]
                product_names = value[2].split('\t')
            else:
                common_name = ''
                product_names = value[1].split('\t')
                
            if common_name.strip() == '':
                common_name = '\N'
            for product_name in product_names:
                if product_name in self.protein_list:
                    self.membership_list.append('%s\t%s' % (product_name, gene_name))
                    gene_set.add((gene_name, common_name, 'gene', self.run_name))
        gene_list = ['\t'.join(list(val)) for val in gene_set]
        self.write_and_close(os.path.join(self.run_dir,
                                          'metacyc_gene.tab_delimited'),
                             gene_list)

    def return_first_value_or_null(self, value):
        value_list = value.split('\t')
        if len(value_list) >= 1:
            if value_list[0].strip() != '':
                return value_list[0]
        return '\N'

    def get_reaction_data(self):
        '''Creates the reaction file and adds associations in the membership_file.'''
        print 'Total proteins:', len(self.protein_list)
        protein_string = self.convert_to_sql_list(self.protein_list)
        # First check find all enzymatic reactions where these are present as enzymes
        # update the reaction list and membership list as well
        self.cursor.execute('select reaction, enzyme from enzrxns where enzyme in %s'
                            % protein_string)
        values = self.cursor.fetchall()
        print 'Obtained by enzrxns:', len(values)
        for value in values:
            reaction_list = value[0].split('\t')
            for reaction in reaction_list:
                self.reaction_list.append(reaction)
                self.membership_list.append('%s\t%s' % (value[1], reaction))
        obtained_values = [value[1] for value in values]
        no_reaction_proteins = list(set(self.protein_list).difference(set(obtained_values)))
        print 'Remaining:', len(no_reaction_proteins)
        # Now check the remaining values in reaction left
        self.cursor.execute('select unique_id, left from reactions')
        values = self.cursor.fetchall()
        left_values = [value[1] for value in values]
        obtained_values = set()
        for protein in no_reaction_proteins:
            for reaction_id, left_val in values:
                if protein in left_val:
                    obtained_values.add(protein)
                    self.reaction_list.append(reaction_id)
                    self.membership_list.append('%s\t%s' % (protein, reaction_id))
        print 'Obtained from reactions left: ', len(obtained_values)
        no_reaction_proteins = list(set(no_reaction_proteins).difference(obtained_values))
        print 'Remaining:', len(no_reaction_proteins)
        # Now check the remaining values in reaction right
        self.cursor.execute('select unique_id, right from reactions')
        values = self.cursor.fetchall()
        right_values = [value[1] for value in values]
        obtained_values = set()
        for protein in no_reaction_proteins:
            for reaction_id, right_val in values:
                if protein in right_val:
                    obtained_values.add(protein)
                    self.reaction_list.append(reaction_id)
                    self.membership_list.append('%s\t%s' % (protein, reaction_id))
        print 'Obtained from reactions right: ', len(obtained_values)
        no_reaction_proteins = list(set(no_reaction_proteins).difference(obtained_values))
        print 'Remaining:', len(no_reaction_proteins)
        # Now check for those proteins which are in complexes that are in a reaction
        # Sometimes the component_of field is just missing. Hence the try...except
        try:
            no_reaction_proteins_sql = self.convert_to_sql_list(no_reaction_proteins)
            self.cursor.execute("select unique_id, component_of from proteins where unique_id in %s and component_of != ''" % no_reaction_proteins_sql)
            values = self.cursor.fetchall()
            component_set = set()
            for value in values:
                for component_of_val in value[1].split('\t'):
                    component_set.add((value[0], component_of_val))
            for protein, complex_name in component_set:
                self.membership_list.append('%s\t%s' % (protein, complex_name))
            component_of_values = [value[1] for value in component_set]
            component_of_values = list(set(component_of_values))
            component_of_values_sql = self.convert_to_sql_list(component_of_values)
            self.cursor.execute("select unique_id, catalyzes from proteins where unique_id in %s and catalyzes !=''"
                                % component_of_values_sql)
            values = self.cursor.fetchall()
            for protein, enz_reaction in values:
                enz_reaction = enz_reaction.split('\t')
                enz_reaction_sql = self.convert_to_sql_list(enz_reaction)
                self.cursor.execute("select unique_id, reaction from enzrxns where unique_id in %s"
                                                % enz_reaction_sql)
                reactions = self.cursor.fetchall()
                for enzrxn_id, reaction_list in reactions:
                    reaction_list = reaction_list.split('\t')
                    for reaction in reaction_list:
                        self.reaction_list.append(reaction)
                        self.membership_list.append('%s\t%s' % (protein, reaction))

            values = [value[0] for value in values]
            obtained_values = [value[0] for value in list(component_set) if value[1] in values]
            print 'Obtained from Complexes: ', len(obtained_values)
            no_reaction_proteins = list(set(no_reaction_proteins).difference(obtained_values))
        except:
            print 'No complex information is stored in the dat files ...'
        print 'Remaining:', len(no_reaction_proteins)
            
    def convert_to_sql_list(self, prot_list):
        protein_string = ["'%s'" % protein for protein in prot_list]
        return "(%s)" % ','.join(protein_string)

    
    def write_reaction_file(self):
        '''Write down the reaction file.'''
        #Remove duplicates
        self.reaction_list = list(set(self.reaction_list))
        print 'Number of reactions in original list: ', len(self.reaction_list)
        reaction_sql = self.convert_to_sql_list(self.reaction_list)
        try:
            self.cursor.execute("""select unique_id, common_name, ec_number, systematic_name
            from reactions where unique_id in %s""" %
                                reaction_sql)
        except:
            #Account for cases without systematic_name
            self.cursor.execute("""select unique_id, common_name, ec_number
            from reactions where unique_id in %s""" %
                                reaction_sql)
        reaction_values = self.cursor.fetchall()
        output = []
        for reaction in reaction_values:
            temp_output = []
            temp_output.append(reaction[0])
            temp_output.append(self.return_first_value_or_null(reaction[1]))
            temp_output.append('reaction')
            temp_output.append(self.run_name)
            temp_output.append(self.return_first_value_or_null(reaction[2]))
            #Account for cases without systematic_name
            if len(reaction) == 4:
                temp_output.append(self.return_first_value_or_null(reaction[3]))
            else:
                temp_output.append('\N')
            output.append('\t'.join(temp_output))
        print 'Number of reactions written to file: ', len(output)
        self.write_and_close(os.path.join(self.run_dir,
                                          'metacyc_reaction.tab_delimited'),
                             output)
            
        
    def get_pathway_data(self, selected_pathway=None):
        '''Creates the pathway file that includes every other pathway or reaction.'''
        if not selected_pathway:
            self.cursor.execute("select unique_id, reaction_list from pathways")
        else:
            self.cursor.execute("select unique_id, reaction_list from pathways where unique_id = '%s'" %
                                selected_pathway)
        values = self.cursor.fetchall()
        for pathway_id, reaction_list in values:
            self.pathway_list.append(pathway_id)
            reactions = reaction_list.split('\t')
            for reaction in reactions:
                self.membership_list.append('%s\t%s' % (reaction, pathway_id))
                if 'RXN' in reaction:
                    self.reaction_list.append(reaction)
                elif 'PWY' in reaction:
                    self.get_pathway_data(reaction)

    def write_pathway_file(self):
        self.pathway_list = list(set(self.pathway_list))
        print 'Number of pathways in original list: ', len(self.pathway_list)
        pathway_sql = self.convert_to_sql_list(self.pathway_list)
        self.cursor.execute("""select unique_id, common_name from pathways where unique_id in %s""" %
                            pathway_sql)
        pathway_values = self.cursor.fetchall()
        output = []
        for pathway in pathway_values:
            temp_output = []
            temp_output.append(pathway[0])
            temp_output.append(self.return_first_value_or_null(pathway[1]))
            temp_output.append('pathway')
            temp_output.append(self.run_name)
            output.append('\t'.join(temp_output))
        print 'Number of pathways written to file: ', len(output)
        self.write_and_close(os.path.join(self.run_dir,
                                          'metacyc_pathway.tab_delimited'),
                             output)

    def write_membership_file(self):
        self.membership_list = list(set(self.membership_list))
        self.write_and_close(os.path.join(self.run_dir,
                                          'metacyc_membership.tab_delimited'),
                             self.membership_list)

    def create_sqlite_database(self, dir_name, database_name, files):
        '''Creates an sqlite_database,
        adds information from the dat files and
        overwrites the database if already present.'''
        print 'Creating temporary sqlite database.'''
        self.create_and_conn_to_sqlite_db(database_name)
        for fileval in files:
            column_name = self.get_column_name(fileval)
            print '%s column in sqlite db started ...' % column_name
            dat_file_dict = self.parse_metacyc_file(fileval)
            self.convert_dict_into_table(dat_file_dict, column_name)
            print '%s column done ...' % column_name
        print 'Sqlite database %s created ...' % database_name

    def write_postgres_sql_file(self):
        '''Creates a new postgres sql file based on the newly created values.'''
        postgres_template = open('/home/awarrier/ohana_repository/bpg/common/parsers/metacyc_parsers/metacyc_creation_template.sql').read()
        cur_dir = {'cur_dir' : os.path.abspath(self.run_dir)}
        output_string = postgres_template % cur_dir
        f = open(os.path.join(self.run_dir,'metacyc_creation.sql'), 'w')
        f.write(output_string)
        f.close()

    def get_accessions(self, input_file):
        '''Get accessions from fasta file.''' 
        headers = [line for line in open(input_file).readlines()
                   if line.startswith('>')]
        accs = [self.get_accession_from_header(header)
                        for header in headers]
        return list(set(accs))

    def get_accession_from_header(self, header):
        '''Get uniprot accession from header regardless of source.'''
        # For fasta files from QuestForOrthologs
        if ':' in header.split()[0]:
            header_start = header.split()[0]
            accession = header_start.split(':')[1]
        elif '|' in header.split()[0]:
            header_start = header.split()[0]
            accession = header_start.split('|')[1]
            if accession.lower() == self.run_name:
                accession = header_start.split('|')[2]
        else:
            accession = header.split()[0]
            #Some headers have a weird space after the > symbol
            if accession.strip() == '>':
                accession = header.split()[1]
        accession = accession.replace('>', '')
        return accession

    def get_sequence(self, blastable_db, accession):
        '''Get sequence from uniprot id'''
        seq = subprocess.Popen(['fastacmd', '-d', blastable_db,
                                '-s', accession],
                               stdout = subprocess.PIPE)
        return seq.stdout.read()

    def get_seguid(self, sequence):
        sequence_handle = StringIO.StringIO(sequence)
        record = SeqIO.parse(sequence_handle, 'fasta').next()
        sequence_handle.close()
        seguid = CheckSum.seguid(record.seq)
        return seguid

    def check_in_uniprot_table_based_on_seguids(self, accession_list):
        '''Returns a dictionary for PF3.0 included proteins in the accession list, where the
        keys are the biocyc protein product accessions. Also writes down a list of unincluded
        proteins.'''
        start_time = time.time()
        unincluded_proteins = []
        included_proteins = {}
        for accession in accession_list:
            sequence = self.get_sequence(self.protein_file, accession)
            seguid = self.get_seguid(sequence)
            accessions = sp.sql_results(cursor,
                                        "select distinct accession from uniprot where seguid=%s;",
                                        (seguid,))
            if len(accessions) == 0:
                unincluded_proteins.append(accession)
            else:
                included_proteins[accession] = []
                for acc in accessions:
                    included_proteins[accession].append(acc[0])
        self.write_and_close(os.path.join(self.run_dir, 'unincluded_proteins'),
                             unincluded_proteins)
        print 'Elapsed: %s' % (time.time() - start_time)
        return included_proteins

    def check_in_uniprot_table_based_on_seguids2(self, accession_list):
        '''Returns a dictionary for PF3.0 included proteins in the accession list, where the
        keys are the biocyc protein product accessions. Also writes down a list of unincluded
        proteins. Replaced previous version to search for seguids in one go. Not much of a time saver
        and more confusing.'''
        start_time = time.time()
        seguid_dict = dict([(accession, self.get_seguid(
            self.get_sequence(self.protein_file, accession)))
                            for accession in accession_list])
        # Get seguids
        seguids = list(set(seguid_dict.values()))
        seguid_sql = self.convert_to_sql_list(seguids)
        # Get uniprot accessions for these seguids in one shot
        seguid_unip_list = sp.sql_results_without_params(cursor,
                                                         """select seguid, accession from
                                                         uniprot where seguid in %s"""
                                                         % seguid_sql)
        # Create a uniprot seguid dictionary
        seguid_unip_dict = {}
        for seguid, unip_acc in seguid_unip_list:
            if not seguid in seguid_unip_dict:
                seguid_unip_dict[seguid] = []
            seguid_unip_dict[seguid].append(unip_acc)
        #Build the included_proteins dictionary
        included_proteins = {}
        for accession, seguid in seguid_dict.items():
            if seguid in seguid_unip_dict:
                included_proteins[accession] = list(set(seguid_unip_dict[seguid]))
        included_proteins_set = set(included_proteins.keys())
        all_proteins = set(accession_list)
        unincluded_proteins = list(all_proteins - included_proteins_set)
        self.write_and_close(os.path.join(self.run_dir, 'unincluded_proteins'),
                             unincluded_proteins)
        print 'Elapsed: %s' % (time.time() - start_time)
        return included_proteins    


    def write_and_close(self, filename, list_vals):
        '''Simple utility to write to a file and close it.'''
        f = open(filename, 'w')
        f.write('\n'.join(list_vals))
        f.write('\n')
        f.close()

    def write_mapping(self):
        '''Writes a file with the id mapping between the metacyc protein and
        uniprot information.'''
        f = open(os.path.join(self.run_dir, 'mapping'), 'w')
        for i in self.uniprot_mapping.items():
            f.write('%s\n' % str(i))
        f.close()
        print 'Mapping file written in the run directory ...'

    def create_mapping(self):
        '''Creates a mapping between uniprot and the protein file within the input
        folder.'''
        accessions = self.get_accessions(self.protein_file)
        print '%s Biocyc protein accessions obtained ...' % len(accessions)
        self.uniprot_mapping = self.check_in_uniprot_table_based_on_seguids(accessions)
        print 'Mapping onto UniProt in database completed ...'
        print 'Matches were found for %s BioCyc proteins ...' % len(self.uniprot_mapping)
        
    
    def main(self):
        '''Runs process for each file in dat_files'''
        # Checks if required files for parsing are present
        files = self.check_if_required_files_present(self.dat_dir)
        if not files:
            print 'Failed for %s ...' % self.run_name
            return
        # Since the membership file is the last one created, do not continue if this file is found.
        membership_file = os.path.join(self.run_dir, 'metacyc_membership.tab_delimited')
        if os.path.exists(membership_file):
            print 'Membership file already exists ...'
            print 'Quitting ...'
            return
        # Create working directory
        self.make_dir(self.run_dir)
        # Creates a mapping string between the uniprot accessions and metacyc_products
        self.create_mapping()
        # Write the mapping string to 'mapping'
        self.write_mapping()
        # Create sqlite database
        self.create_sqlite_database(self.dat_dir,
                                    os.path.join(self.run_dir,
                                                 '%s.db' % self.run_name),
                                    files)
        # Create protein file
        self.create_protein_file()
        # Create gene file
        self.create_gene_file()
        # Get reaction and pathway data
        self.get_reaction_data()
        self.get_pathway_data()
        # Write this data and the membership file
        self.write_reaction_file()
        self.write_pathway_file()
        self.write_membership_file()
        # Close connection to sqlite database
        self.close_connection()
        # Write postgres sql file
        self.write_postgres_sql_file()
        
if __name__ == "__main__":
    dat_dir = sys.argv[1]
    start_time = time.time()
    mp = MetacycParser(dat_dir)
    mp.main()
    #mp.write_postgres_sql_file()
    print 'Time taken: %s seconds' % (time.time() - start_time)
