#!/usr/bin/env python2.7
'''Simple utility that goes through all the folders in the data directory
created by metacyc.py and collects the different tab delimited files into
a single one into the collection directory.'''

import os
import sys

class MetacycCollector:
    def __init__(self, version):
        '''No input arguments at this point except the version number of biocyc.
        '''
        self.version = version
        self.protein_set = set()
        self.gene_set = set()
        self.reaction_set = set()
        self.pathway_set = set()
        self.membership_set = set()
        self.output_dir = 'collection'
        try:
            os.mkdir(self.output_dir)
        except:
            pass

    def check_if_data_directory_present(self):
        '''Run in the previous working directory, checks for the presence of the
        data directory and reports the number of directories within.
        '''
        if not os.path.exists('data'):
            print 'Data directory is absent ...'
            print 'Exiting ...'
            sys.exit()
        else:
            num_dirs = len(os.listdir('data'))
            print 'There are %s directories or databases in the data dir ...' % num_dirs

    def collect_files(self, ind_filename):
        '''Collects all protein files using a simple cat command.'''
        os.system('for file in data/*/%s; do cat $file >> collection/%s; echo "" >> collection/%s;done' %
                  (ind_filename, ind_filename, ind_filename))

    def clean_protein_file(self):
        '''Cleans the protein file to remove any escape characters.'''
        print 'Cleaning protein file ...'
        f = open('collection/metacyc_protein.tab_delimited').readlines()
        f2 = open('collection/metacyc_protein.tab_delimited', 'w')
        for line in f:
            if not line.strip() == '':
                line_split = line.split('\t')
                if not line_split[1] == '\N':
                    line_split[1] = line_split[1].replace('\\', '')
                    line_split[1] = line_split[1].replace("'", "")
                f2.write('\t'.join(line_split))
        f2.write('\n')
        f2.close()

    def clean_reaction_file(self):
        '''Cleans the reaction file to remove any escape characters.'''
        print 'Cleaning reaction file ...'
        f = open('collection/metacyc_reaction.tab_delimited').readlines()
        f2 = open('collection/metacyc_reaction.tab_delimited', 'w')
        for line in f:
            if not line.strip() == '':
                line_split = line.split('\t')
                line_split[1] = line_split[1].decode("ascii", "ignore")
                line_split[5] = line_split[5].decode("ascii", "ignore")
                f2.write('\t'.join(line_split))
        f2.write('\n')
        f2.close()

    def remove_newlines(self):
        '''Remove any newlines.'''
        for filename in ['gene', 'protein', 'reaction', 'pathway', 'membership']:
            print 'Removing newlines in %s ...' % filename
            f = open('collection/metacyc_%s.tab_delimited' % filename).readlines()
            f2 = [line for line in f if not line.strip() == '']
            f = open('collection/metacyc_%s.tab_delimited' % filename, 'w')
            f.write(''.join(f2))
            f.close()
        
    def clean_membership_file(self):
        '''Removes duplicates from the membership file.'''
        print 'Cleaning membership file ...'
        f = open('collection/metacyc_membership.tab_delimited').readlines()
        f_set = set(f)
        f2 = open('collection/metacyc_membership.tab_delimited', 'w')
        f2.write(''.join(list(f_set)))
        f2.close()
        
        
    def list_filecount(self, ind_filename):
        '''Just a count of the members in this file.'''
        print os.system('wc -l %s' % ind_filename)

    def write_postgres_sql_file(self):
        '''Creates a new postgres sql file based on the newly created values.'''
        postgres_template = open('/home/awarrier/ohana_repository/bpg/common/parsers/metacyc_parsers/metacyc_creation_template.sql').read()
        cur_dir = {'cur_dir' : os.path.abspath(self.output_dir)}
        output_string = postgres_template % cur_dir
        f = open(os.path.join(self.output_dir,'metacyc_creation.sql'), 'w')
        f.write(output_string)
        f.close()
        
    def main(self):
        '''Main function.'''
        self.check_if_data_directory_present()
        # Collect protein files
        print 'Collecting protein files ...'
        self.collect_files('metacyc_protein.tab_delimited')
        # Collect gene files
        print 'Collecting gene files ...'
        self.collect_files('metacyc_gene.tab_delimited')

        # Collect reaction files
        print 'Collecting reaction files ...'
        self.collect_files('metacyc_reaction.tab_delimited')

        # Collect pathway files
        print 'Collecting pathway files ...'
        self.collect_files('metacyc_pathway.tab_delimited')

        # Collect membership files
        print 'Collecting membership files ...'
        self.collect_files('metacyc_membership.tab_delimited')
        # Clean the files
        self.clean_protein_file()
        self.clean_reaction_file()
        self.clean_membership_file()
        self.remove_newlines()
        # write version file
        f = open('collection/metacyc_version.tab_delimited', 'w')
        f.write(self.version)
        f.close()
        #write sql file
        self.write_postgres_sql_file()
        
        
        
if __name__ == "__main__":
    version = sys.argv[1]
    mc = MetacycCollector(version)
    mc.main()

                

    
