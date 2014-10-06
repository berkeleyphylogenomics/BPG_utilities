#!/usr/bin/env python
'''This script (actually a pipeline) takes a fasta file with protein information, obtains best books (both PFAM and GHG) for each,
creates download files (HMMs, MSAs, trees, summaries)
moves them to the download location,
writes a config file that is then displayed on the website.
Please read the main function to understand the workflow of this pipeline

Pre-requisites:
Multiple changes have occured in this script.
Previous versions required a PFAM-A hmm libarary. 

WARNING:
DO NOT RUN THIS SCRIPT AS A SHARD JOB. The script itself spawns shard jobs several times.
If the script is non-responsive for too long, check "qstat" to find out if the shards are
still running. If they have completed and there are no errors, press Ctrl-C to continue.
The script will resume and assume that you are satisfied with the shard job run.
'''
import os
import sys
import subprocess
import re
import pwd
import shlex
import shutil
import logging
import time

# Setting the logger so that information is written onto files
# This code is from
# http://www.mechanicalcat.net/richard/log/Python/Simple_usage_of_Python_s_logging_module
logging.basicConfig()


class GenomePageCreator:
    '''Main class called with an input config dictionary.'''
    def __init__(self, config_dict):
        self.start_time = time.time()
        #First read from config dictionary and populate instance variables.
        #Would have been simpler to just use the config dictionary as an instance variable
        #itself. However, I'm laying it out this way so that the individual elements
        #are easy to follow.
        self.genome_name = config_dict['genome_name']
        self.run_name = config_dict['genome_short_name']
        self.fasta_file = config_dict['source_fasta_file']
        self.blastable_database_directory = config_dict['blastable_database_directory']
        self.blastable_database = config_dict['blastable_database']
        self.download_location = config_dict['download_location']
        self.out_dir = config_dict['output_dir']
        self.template_file = config_dict['template_for_configuration']
        self.config_file_location = config_dict['configuration_file_location']
        self.num_shards = config_dict['number_of_shards']
        self.error_ignore = config_dict['ignore_errors']
        self.ohana_repository = config_dict['your_ohana_repository_location']
        # The below values are directly carried over to the output
        self.output_config_dictionary = {}
        self.output_config_dictionary['genome_name'] = config_dict['genome_name']
        self.output_config_dictionary['source_reference'] = config_dict['source_reference']
        self.output_config_dictionary['source_link'] = config_dict['source_link']
        # Logging, uniprot lists and some empty lists.
        self.setup_logger()
        self.log(self.run_name)
        self.uniprot_list = self.write_uniprot_list_from_fasta_file()
        self.pfam_uniprot_list = self.write_uniprot_list_for_pfam_shards()
        self.pfam_log = 'PFAM_out'
        self.final_downloads = []
        self.write_config_file_to_config_directory = True
        
    def setup_logger(self):
        '''Set logging so that replacing print statements with self.log
        writes logs to a file.
        This general purpose function can be moved to a common location and used from there.
        '''
        self.logger = logging.getLogger('create_downloads')
        hdlr = logging.FileHandler('%s.log' % self.run_name)
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        self.logger.addHandler(hdlr)
        self.logger.setLevel(logging.INFO)

    def activate_environment(self):
        '''Activates the old environment.'''
        os.system('source /clusterfs/ohana/software/ohana_environment.sh')
        os.system('development_old')
        os.system('development')
        os.system('development_old')
        
    def log(self, string):
        '''Actual print statement replacement.'''
        self.logger.info(string)

    def get_uniprot_accession_from_header(self, header):
        '''Get uniprot accession from header regardless of source.'''
        # For fasta files from QuestForOrthologs
        if ':' in header.split()[0]:
            header_start = header.split()[0]
            uniprot_accession = header_start.split(':')[1]
        elif '|' in header.split()[0]:
            uniprot_accession = header.split('|')[1]
        return uniprot_accession

    def write_uniprot_list_from_fasta_file(self):
        '''Get uniprot_ids from fasta file.''' 
        headers = [line for line in open(self.fasta_file).readlines()
                   if line.startswith('>')]
        uniprot_accs = [self.get_uniprot_accession_from_header(header)
                        for header in headers]
        uniprot_file = '%s_uniprot_list' % self.run_name
        self.write_to_file_and_close(uniprot_file,
                                     '\n'.join(uniprot_accs))
        self.log('\tUniProt list has been written to %s.' % uniprot_file)
        return uniprot_file
        
    def write_uniprot_list_for_pfam_shards(self):
        '''Write a file that includes the uniprot list and locations
        of blastable_db and pfam_hmm_scan_outs. Later changed to include blastable databases
        only'''
        orig_list = open(self.uniprot_list).readlines()
        orig_list = [line.strip() for line in orig_list]
        self.output_config_dictionary['num_genes'] = len(orig_list)
        if len(orig_list) < 5000:
            self.num_shards = 300
        new_list = ['%s\t%s' % (line, self.blastable_database)
                    for line in orig_list]
        self.write_to_file_and_close('pfam_uniprot_list', '\n'.join(new_list))
        self.log('\tUniProt list for Pfam shards has been written.')
        return 'pfam_uniprot_list'
    
    def run_gzip_on_files_in_folder(self):
        '''Run gzip on files in download folder'''
        for file_name in self.final_downloads:
            self.log('\tStart gzipping %s' % file_name)
            cmd = 'gzip %s' % os.path.join('final_downloads',
                                            file_name)
            self.run_subprocess_on_cmd(cmd)

    def get_filesize(self, filename):
        '''Getting human readable filesize'''
        num = os.path.getsize(filename)
        for x in ['bytes','KB','MB','GB','TB']:
            if num < 1024.0:
                return "%3.2f%s" % (num, x)
            num /= 1024.0

    def move_files_to_download_location(self):
        '''If the download location already exists, dont overwrite.
        Just print a message.
        '''
        final_location_for_downloads = os.path.join(self.download_location,
                                                    self.run_name)
        if os.path.exists(final_location_for_downloads):
            self.log('\tFinal downloads location %s already exists. \
            \n\tNot moving files from the downloads folder')
            self.write_config_file_to_config_directory = False
        else:
            os.mkdir(final_location_for_downloads)
            for filename in self.final_downloads:
                gz_file = '%s.gz' % filename
                orig_location = os.path.join('final_downloads',
                                             gz_file)
                final_loc = os.path.join(final_location_for_downloads,
                                         gz_file)
                shutil.copyfile(orig_location, final_loc)
            #Also move the source fasta file to this location
            source_filename = os.path.basename(self.fasta_file)
            final_loc = os.path.join(final_location_for_downloads,
                                     source_filename)
            shutil.copyfile(self.fasta_file, final_loc)
            self.output_config_dictionary['source_fasta_file'] = final_loc
                
    def update_dictionary_for_config(self):
        '''Update config_dictionary to include more info. Basically the filesize and the stats.'''
        for download_type in ['alignments', 'ml_trees',
                              'nj_trees', 'summaries', 'hmms']:
            for book_type in ['bestGHG', 'bestPFAM']:
                download_file = '%s/%s/%s_%s.%s.gz' % (self.download_location,
                                                         self.run_name, self.run_name,
                                                         book_type, download_type)
                book_index = '%s_%s' % (book_type, download_type)
                self.output_config_dictionary['%s_size' % book_index] = self.get_filesize(
                    download_file)
        self.get_GHG_stats()
        self.get_PFAM_stats()
        self.get_all_stats()
                
    def get_GHG_stats(self):
        ghg_list = open('%s_GHG_out' % self.run_name).readlines()
        protein_list = [x for x in ghg_list if not (x.endswith(',')
                                                    or x.endswith(',0'))]
        protein_list = [x for x in protein_list if 'bpg' in x]
        unique_proteins = [x.split(',')[0] for x in protein_list[1:]]
        self.unique_GHG_proteins = set(unique_proteins)
        f = open('%s_GHG_proteins' % self.run_name, 'w')
        f.write('\n'.join(list(self.unique_GHG_proteins)))
        f.close()
        num_genes_in_GHG_pct = (
            float(len(self.unique_GHG_proteins)*100)/self.output_config_dictionary['num_genes'])
        num_genes_in_GHG_pct = "%0.2f" % num_genes_in_GHG_pct
        self.output_config_dictionary['num_genes_in_GHG_pct'] = num_genes_in_GHG_pct
        self.output_config_dictionary['num_genes_in_GHG'] = len(self.unique_GHG_proteins)
        self.log('No of proteins in GHGs = %s' % len(self.unique_GHG_proteins))
        self.log('pct in GHGs = %s' % num_genes_in_GHG_pct)

    def get_PFAM_stats(self):
        '''Gets pfam book stats from the output files
        '''
        pfam_list = open('%s_PFAM_out' % self.run_name).readlines()
        protein_list = [x for x in pfam_list if x.startswith('Best Book')]
        unique_proteins = [x.split('\t')[1] for x in protein_list]
        self.unique_pfam_proteins = set(unique_proteins)
        f = open('%s_PFAM_proteins' % self.run_name, 'w')
        f.write('\n'.join(list(self.unique_pfam_proteins)))
        f.close()
        num_genes_in_pfam_pct = (
            float(len(self.unique_pfam_proteins)*100)/self.output_config_dictionary['num_genes'])
        num_genes_in_pfam_pct = "%0.2f" % num_genes_in_pfam_pct
        self.output_config_dictionary['num_genes_in_PFAM_pct'] = num_genes_in_pfam_pct
        self.output_config_dictionary['num_genes_in_PFAM'] = len(self.unique_pfam_proteins)
        self.log('No of proteins in PFAMs = %s' % len(self.unique_pfam_proteins))
        self.log('pct in PFAMs = %s' % num_genes_in_pfam_pct)

    def get_all_stats(self):
        '''Get stats for all genes.'''
        no_in_any = len(self.unique_pfam_proteins.union(
                                self.unique_GHG_proteins))
        self.log('No of proteins convered in any book = %s' % no_in_any)
        num_genes_in_any_pct = float(no_in_any * 100)/self.output_config_dictionary['num_genes']
        num_genes_in_any_pct = "%0.2f" % num_genes_in_any_pct        
        self.log('pct in any = %s' % num_genes_in_any_pct)
        self.output_config_dictionary['num_genes_in_any_pct'] = num_genes_in_any_pct
        self.output_config_dictionary['num_genes_in_any'] = no_in_any

    def write_to_file_and_close(self, filename, write_string):
        '''Write to a file and close it.'''
        f = open(filename, 'w')
        f.write(write_string)
        f.close()
    
    def create_config_file(self):
        '''Create config file and save within final_downloads and in the django
        accessible location.'''
        temp_string = open(self.template_file).read()
        final_string = temp_string % self.output_config_dictionary
        final_string = final_string.replace(self.download_location,
                                            'http://makana-static.berkeley.edu/genomes')
        self.write_to_file_and_close('final_downloads/config.txt', final_string)
        genome_config_dir = os.path.join(self.config_file_location, self.run_name)
        self.recreate_directory(genome_config_dir)
        self.write_to_file_and_close(os.path.join(genome_config_dir,
                                                  'config.txt'),
                                     final_string)
    
    def cleanup(self):
        '''Remove any intermediate files.'''
        for folder in ['intermediates', 'outputs', 'scripts']:
            try:
                shutil.rmtree(folder)
            except:
                pass
        #Remove shard logs
        for filename in os.listdir('.'):
            if filename.startswith('write_'):
                os.remove(filename)

    def should_blastdb_be_done(self):
        '''Check if blastable database is present
        '''
        if not os.path.exists(self.blastable_database + '.pin'):
            self.log('\tBlastable database for %s is absent.' % self.run_name)
            return True
        else:
            self.log('\tBlastable database for %s found.' % self.run_name)
            return False

    def create_new_fasta_file(self):
        '''Avoid issues with fastacmd not working with the created blastable db.
        Just replacing the entire start chars (between > and the colon or pipe)
        with >lcl|'''
        file_list = open(self.fasta_file).readlines()
        self.lcl_fasta_file = 'lcl_fasta_file_for_%s' % self.run_name
        final_list = []
        for line in file_list:
            if line.startswith('>'):
                if ':' in line.split()[0]:
                    line = re.sub('^.*?:', '>lcl|', line)
                elif '|' in line.split()[0]:
                    line = re.sub('^.*?\|', '>lcl|', line)
            final_list.append(line)
        self.write_to_file_and_close(self.lcl_fasta_file, ''.join(final_list))
        self.log('\tlcl replaced file written.')
        
    def create_blastable_db(self):
        '''Create blastable database with input fasta file.'''
        self.create_new_fasta_file()
        cmd = 'formatdb -i %s -n %s -o T' % (self.lcl_fasta_file,
                                             self.run_name)
        self.run_subprocess_on_cmd(cmd)
        self.log('\tCreated blastable database using formatdb')
        for b_file in os.listdir('.'):
            if b_file.startswith(self.run_name):
                if not (b_file.endswith('.log') or
                        'uniprot' in b_file):
                    shutil.move(b_file, self.blastable_database_directory)
        self.log('\tBlastable files have been move to %s' % self.blastable_database_directory)
        
    def run_hmm_scan_on_pfamA(self):
        '''Running pfam hmm scan on the fasta file. This was done in a previous version of
        the script and left here only for legacy reasons.
        '''
        cmd = 'hmm3scan --cut_ga %s %s' % (PFAM_HMM_LIBRARY, self.fasta_file)
        run_output = self.run_subprocess_on_cmd(cmd)
        self.write_to_file_and_close(PFAM_HMM_SCAN_OUT % self.run_name, run_output)
        self.log('\tHmm scan done and output moved to %s' % (PFAM_HMM_SCAN_OUT
                                                          % self.run_name))


    def recreate_directory(self, folder_name):
        '''Remove a directory and then create a new one again.'''
        try:
            shutil.rmtree(folder_name)
        except:
            pass
        os.mkdir(folder_name)
        
    def make_folders_if_none(self):
        folders = ['outputs', 'scripts', 'outputs_best_GHG',
                   'downloads', 'final_downloads', 'intermediates']
        for folder in folders:
            self.recreate_directory(folder)
        
    def check_if_shard_is_completed(self, directory, num_files):
        try:
            while (len(os.listdir(directory)) < int(num_files)):
                pass
        except KeyboardInterrupt:
            return

    def exit_if_shard_failed(self, filename):
        temp_file = open(filename,'r').read()
        if 'Failed!' in temp_file: # check if there is failed job
            self.log('\tError in %s' % filename)
            if not self.error_ignore:
                sys.exit()
        
    def run_subprocess_on_cmd(self, cmd):
         '''Run subprocess on a command, instead of os.system'''
         args = shlex.split(cmd)
         job = subprocess.Popen(args,stdout=subprocess.PIPE)
         return job.stdout.read()

    def run_cluster_shards(self, input_file, script):
        '''Run cluster_shards on the selected script.'''
        cmd = 'cluster_shards_library -i %s -c %s -s %s -n get_bestPFAMBook -d /home/yshen/python/right_path' % (input_file, script, self.num_shards)
        self.log('Shard command: %s' % cmd)
        self.run_subprocess_on_cmd(cmd)
        
    def get_best_PFAM_book(self):
        '''Get the best PFAM books'''
        #Empty the outputs directory, run pfam best book gathering on the cluster
        #concatenate the outputs and get bpg accessions of best books.
        output_file_name = '%s_bestPFAMBooks.bpgs' % self.run_name
        self.recreate_directory('outputs')
        run_script = ('%s/bpg/common/utils/get_best_Pfam_byAccession.py'
                      % self.ohana_repository)
        self.run_cluster_shards(self.pfam_uniprot_list, run_script)
        self.log('\tStarting pfam best book collection. (Shard job)')
        self.check_if_shard_is_completed('outputs', self.num_shards)
        self.log("\tShard job finished\n\tProcessing the ouputs...")
        pfam_out_file = '%s_PFAM_out' % self.run_name
        os.system('cat outputs/*.out > %s' % pfam_out_file)
        self.exit_if_shard_failed(pfam_out_file)
        self.log('\tNo errors found in the shard job.')
        pfam_best_book_lines = [line for line in open(pfam_out_file).readlines()
                                if line.startswith('Best Book:')]
        pfam_best_book_lines = [line.split('\t') for line in pfam_best_book_lines]
        pfam_best_families = []
        for each_best_book in pfam_best_book_lines:
            pfam_best_families += [item for item in each_best_book if item.startswith('bpg')]
        pfam_best_families = set(pfam_best_families)
        self.log('\tThere are %s pfam books.' %
                           len(list(pfam_best_families)))
        self.output_config_dictionary['num_pfam_books'] = len(list(pfam_best_families))
        self.write_to_file_and_close(output_file_name, '\n'.join(pfam_best_families))
        self.log('\tPFAM best books obtained and written.')
        return output_file_name

    def get_best_GHG_book(self):
        '''Get the best GHG books'''
        #remove_outputs_folder_if_exists()
        self.recreate_directory('outputs')
        run_script = ('%s/bpg/common/utils/get_best_GHG_byAccession.py'
                              % self.ohana_repository)
        self.run_cluster_shards(self.uniprot_list, run_script)
        self.log('\tStarting GHG best book collection. (Shard job)')
        self.check_if_shard_is_completed('outputs', self.num_shards)
        os.system('cat outputs/*.out > %s_GHG_out' % self.run_name)
        # put individual output files into one file, and extract the bpg accessions
        self.log("\tShard job finished\n\tProcessing the ouputs...")
        self.exit_if_shard_failed('%s_GHG_out' % self.run_name)
        self.log('\tNo errors found in the shard job.')
        bpg_list_file = '%s_bestGHGBooks.bpgs' % self.run_name
        GHG_out = open('%s_bestGHGBooks.out' % self.run_name,'w')
        bpgs = []
        for each_out in os.listdir('outputs'):
            filename = os.path.join('outputs', each_out)
            self.exit_if_shard_failed(filename)
            for line in open(os.path.join('outputs', each_out),'r').readlines():
                if 'bpg' in line:
                    GHG_out.write('%s' % line)
                    if re.search('(bpg\d{7})',line):
                        bpg_acc = re.search('(bpg\d{7})',line)
                        bpgs.append(str(bpg_acc.group(1)))
        bpg_nr = set(bpgs)
        self.write_to_file_and_close(bpg_list_file, '\n'.join(list(bpg_nr)))
        self.output_config_dictionary['num_GHG_books'] = len([line for line in
                                                       open(bpg_list_file).readlines()
                                                       if not line.strip() == ''])
        GHG_out.close()
        GHG_dir = "%s_best_GHG" % self.out_dir
        try:
            shutil.move(self.out_dir, "%s_best_GHG" % self.out_dir)
        except:
            pass
        self.log('\tThere are %s GHG books.' %
                           len(list(bpg_nr)))
        self.log('\tGHG best books obtained and written.')
        return(bpg_list_file)

    def print_number_of_bpg_accessions(self, filename):
        '''Just run a basic grep and get the number of bpg accession lines
        in filename'''
        cmd = "grep '^bpg' %s" % filename
        num = self.run_subprocess_on_cmd(cmd)
        num = len([line for line in num.split()
                   if not line.strip() == ''])
        self.log('\t%s : %s bpg lines' % (filename,
                                                  num))

    def concatenate_downloads(self, type):
        '''Concatenate individual files to specific downloads.'''
        final_downloads = []
        for download_type in ['alignments', 'ml_trees',
                              'nj_trees', 'summaries']:
            download_name = '%s%s_%s' % (type, self.run_name, download_type)
            final_download_name = '%s_%s.%s' % (self.run_name, type, download_type)
            #concatenate using 'cat', use grep to get the number of bpgs
            # and then move to final_downloads folder, after deleting the
            # individual files
            self.log('\tConcatenating %s' % final_download_name)
            os.system('cat downloads/%s_* > downloads/%s' % (
                download_name, final_download_name))
            self.log('\tRemoving intermediate download files.')
            os.system('rm downloads/%s_*' % download_name)
            self.print_number_of_bpg_accessions('downloads/%s' % final_download_name)
            self.log('\tMoving %s %s file to final_downloads folder.' % (type, download_type))
            shutil.move('downloads/%s' % final_download_name,
                        'final_downloads/%s' % final_download_name)
            final_downloads.append(final_download_name)
            self.output_config_dictionary['%s_%s' % (type, download_type)] = ('%s/%s/%s.gz' %
                                                                       (self.download_location,
                                                                        self.run_name,
                                                                        final_download_name))
        #Separately for hmms since no concatenation is required
        shutil.move('downloads/%s%s.hmms' % (type, self.run_name),
                        'final_downloads/%s_%s.hmms' % (self.run_name, type))
        final_download_name = '%s_%s.%s' % (self.run_name, type, 'hmms')
        self.output_config_dictionary['%s_%s' % (type, 'hmms')] = ('%s/%s/%s.gz' %
                                                                       (self.download_location,
                                                                        self.run_name,
                                                                        final_download_name))
        final_downloads.append(final_download_name)
        self.final_downloads += final_downloads
        self.log('\tFinished concatenating downloads')
        
    def create_download(self, bpg_list, type):
        '''Create all downloads given the bpg_list and type'''
        tag=type+self.run_name
        for method in ['ml','nj']:
            cmd='shard_job.py write_trees_to_filesystem \
            "%s/bpg/common/utils/write_trees_to_shard.py \
            -d %s -m %s %s" 30 -lwalltime=72:00:00 -q library' % (
            self.ohana_repository, tag, method, bpg_list)
            self.run_subprocess_on_cmd(cmd)

        # Create the download files for alignments, summaries
        for script in ['alignments','summaries']:
            script_name=('%s/bpg/common/utils/write_%s_to_shard.py'
                         % (self.ohana_repository, script))
            cmd = 'shard_job.py write_%s_to_filesystem "%s -d %s %s" 30 -lwalltime=72:00:00 -q library' % (
                script, script_name, tag, bpg_list)
            self.run_subprocess_on_cmd(cmd)
        # Create HMM downloads
        script_name='%s/bpg/common/utils/write_hmms_to_shard.py' % self.ohana_repository
        cmd = '%s -d %s -i %s' % (script_name, tag, bpg_list)
        self.run_subprocess_on_cmd(cmd)
        #Check if download shard job is done. This should generate 30 files for
        #nj, ml, alignments, and summaries and one file for hmm.
        self.log('\tShard jobs running.')
        self.check_if_shard_is_completed('downloads', 121)
        self.log("\tFinish creating download files for %s." % type)
        self.concatenate_downloads(type)
        
    def main(self):
        '''Main function.'''
        self.activate_environment()
        self.log(sys.path)
        if self.should_blastdb_be_done():
            self.create_blastable_db()
        self.log('Creating directories.')
        self.make_folders_if_none()
        self.log("Start search for best PFAM books")
        pfam_list = self.get_best_PFAM_book()
        self.log("Start creating download files for best PFAM books")
        self.create_download(pfam_list, 'bestPFAM') 
        self.log("Start search for best GHG books")
        ghg_list = self.get_best_GHG_book()
        self.log("Start creating download files for best GHG books")
        self.create_download(ghg_list, 'bestGHG')
        self.log('Running gzip on files.')
        self.run_gzip_on_files_in_folder()
        self.log('Moving gz files.')
        self.move_files_to_download_location()
        if not self.write_config_file_to_config_directory:
            self.log('Exiting. \
            \n A new run_name is needed to continue \
            \n Configuration file is not automatically created. \
            \n Use the template to create this')
            sys.exit()
        self.log('Updating dictionary')
        self.update_dictionary_for_config()
        self.log('Creating a config file')
        self.create_config_file()
        self.log('Cleaning up. Time taken so far: %s' % (float(time.time() - self.start_time)/60))
        self.cleanup()
#Class ends here

def read_config_file(filename):
    '''Read config file and get elements.'''
    filename = os.path.abspath(filename)
    config_dict = open(filename).read()
    config_dict = str(config_dict.strip())
    config_dict = eval(config_dict)
    return config_dict
    

if __name__=='__main__':
    if len(sys.argv)<2:
        print len(sys.argv)
        print """
        Usage %s <Path_to_run_config_file>

        Run this code in an empty directory or in one that you have write access to. A number of intermediate files are created.
        
        Sample contents of the run config file for Eschericia coli.

        {
        'source_fasta_file' : '/clusterfs/ohana/external/genomes/QuestForOrthologs/current/83333_escherichia_coli.fasta',
        'source_reference' : 'Quest For Orthologs, Release 2011_04',
        'source_link' : 'http://www.ebi.ac.uk/reference_proteomes/',
        'output_dir' : 'outputs',
        'number_of_shards' : '150',
        'genome_name' : 'Escherichia coli K-12',
        'genome_short_name' : 'e_coli',
        'your_ohana_repository_location' : '/home/awarrier/ohana_repository',
        'ignore_errors' : 'yes',
        'blastable_database_directory' : '/clusterfs/ohana/external/blastable_db_for_genome_pages/',
        'blastable_database' : '/clusterfs/ohana_repository/external/blastable_db_for_genome_pages/e_coli',
        'download_location' : '/clusterfs/ohana/software/webserver/bigfile_static/genomes',
        'template_for_configuration' : '/clusterfs/ohana/bpg/pfacts/genomes/config.template',
        'configuration_file_location' : '/clusterfs/ohana/bpg/pfacts/genomes/'
        }
        """ % sys.argv[0]
        sys.exit(0)
    start_time = time.time()
    config_dict = read_config_file(sys.argv[1])
    gpc = GenomePageCreator(config_dict)
    gpc.main()
    "Elapsed time: %s min" % (float(time.time() - start_time)/60)
