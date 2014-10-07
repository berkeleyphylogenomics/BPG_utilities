'''
Functions for dealing with BLAST.
'''
import subprocess, tempfile

BLASTDB = '/clusterfs/ohana/external/UniProt/2014-07-31/blastdb/prot'
TEMP_DIR = '/clusterfs/vasudha/software/webserver/temp/'

def retrieve_fasta_for_accession_list(accession_list, blastdb=BLASTDB, temp_dir=TEMP_DIR):
    fasta = ''
    accession_string = '\n'.join(accession_list)
    acc_file = tempfile.NamedTemporaryFile(dir=temp_dir)
    acc_file.write(accession_string)
    args = ['/clusterfs/ohana/software/bin/blastdbcmd', '-db', 
            blastdb, '-dbtype', 'prot', '-entry_batch', acc_file.name]
    acc_file.seek(0)
    process = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    fasta = process.communicate()[0]
    acc_file.close()
    return fasta
