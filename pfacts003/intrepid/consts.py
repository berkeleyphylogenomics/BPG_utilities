##################################################################
# TIMEOUT INTERVALS
##################################################################
PROGRESS_PAGE_TIMEOUT_INTERVAL = 10000
MAIN_PAGE_TIMEOUT_INTERVAL = 10000
##################################################################
# STATE MACHINE STUFF
##################################################################
# This section has the status id part.  In order to allow easy insertion of 
# pipeline stages, we will allow the status strings to be constants defined here
# The integer corresponds to the key that is in the database entry for that stage

# To add a new stage:
# You should really stop the daemon and make sure no jobs are executing before doing this.
#
# 1. Create the stage's row in the database, make sure to increase the total_pipeline_stages accordingly
# 2. For all of the other stage rows (ids of these below), increment the total_pipeline_stages and for
#    stages after the inserted stage, increment their pipeline_stage also
# 3. Add the global name of the stage below and set it equal to the id of the row you created in 1.

# Stage 0 - Job was created.
JOB_CREATED = 1
# Stage 1 - Job can be picked up by the daemon now
WAITING_FOR_SUBMISSION = 2
# Stage 2 - Job was picked up by the daemon, now it's in the queue
JOB_SUBMITTED = 3
# PIPELINE STAGES
# The control is now in the intrepid_pipeline.py file
SCAN_PFAMS = 4
# Stage 3 - First gather homologs
GATHER_HOMOLOGS = 5
# Stage 4 - Masking msa of query and homologs
PROCESS_MSA = 6
# Stage 5 - Building tree of query and homologs
BUILD_TREE = 7
# Stage 6 - Running INTREPID on tree and MSA
RUN_INTREPID = 8
# Stage 7 - Running DISCERN
RUN_DISCERN = 10
# Stage 7 - Wrap up.
JOB_DONE = 9
##################################################################
# DOMAIN LENGTH CONSTRAINTS
##################################################################
MINIMUM_PFAM_DOMAIN_LENGTH = 50 # aa
##################################################################
# PDB CONSTANTS
##################################################################
GENERIC_PDB_ID = 'BPG0'
##################################################################
# THROTTLING
##################################################################
ANONYMOUS_JOB_LIMIT = 10
AUTHENTICATED_USER_THROTTLE_MESSAGE = "You have submitted %d Intrepid jobs since midnight PST.  You must wait until tomorrow before submitting again.  Please email phylofacts.webmaster@gmail.com to increase your daily job submission quota."
ANONYMOUS_USER_THROTTLE_MESSAGE = 'We have received %d Intrepid job submissions from this IP address in the last 24 hours.  This is the maximum allowed for non-registered users.  To become a registered user and increase your quota, please click <a href="/phylofacts/user_create/" target="_blank">here</a>.  Otherwise, you can submit new jobs after midnight PST.'
##################################################################
# FASTA VALIDATION
##################################################################
ACCEPTED_SEQUENCE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
NO_FASTA_MESSAGE = 'FASTA input required'
NO_PDB_MESSAGE = 'PDB input required'
##################################################################
# SEQUENCE LENGTH VALIDATION
##################################################################
MINIMUM_INPUT_LENGTH = 50
MAXIMUM_INPUT_LENGTH = 2000
INPUT_SEQUENCE_TOO_SHORT_MESSAGE = 'Sequence too short.  The PhyloFacts Intrepid Webserver does not accept inputs &lt;%d amino acids' % MINIMUM_INPUT_LENGTH
INPUT_SEQUENCE_TOO_LONG_MESSAGE = 'Sequence too long.  The PhyloFacts Intrepid Webserver does not accept input &gt;%d amino acids' % MAXIMUM_INPUT_LENGTH
##################################################################
# EMAIL VALIDATION
##################################################################
EMAIL_INVALID_MESSAGE = 'Invalid email'
##################################################################
# EMAIL SUBJECT VALIDATION
##################################################################
MAXIMUM_EMAIL_SUBJECT_CHARACTERS = 100
EMAIL_SUBJECT_TOO_LONG_MESSAGE = 'Email subject must be &lt;%d characters' % MAXIMUM_EMAIL_SUBJECT_CHARACTERS
##################################################################
# GAPPY COLUMN CONSTRAINTS
##################################################################
GAPPY_COLUMN_CONDITION = 0.5 # if greater than 50% gaps, remove from alignment before tree building
#GAPPY_COLUMN_CONDITION = -1 # remove a column with any gaps in the seed
##################################################################
# MINIMUM PID TO QUERY FOR QUITTING
##################################################################
MINIMUM_PID_TO_QUERY = 0.15
##################################################################
# INTREPID INCLUDE DIRECTORY
##################################################################
INTREPID_DATA_DIR = '/clusterfs/ohana/software/intrepid/data/'
##################################################################
# AMINO ACID DICTIONARY
##################################################################
AA_DICT = {
    'A': {'3':'Ala','full':'Alanine'},
    'R': {'3':'Arg','full':'Arginine'},
    'N': {'3':'Asn','full':'Asparagine'},
    'D': {'3':'Asp','full':'Aspartic acid'},
    'C': {'3':'Cys','full':'Cysteine'},
    'E': {'3':'Glu','full':'Glutamic acid'},
    'Q': {'3':'Gln','full':'Glutamine'},
    'G': {'3':'Gly','full':'Glycine'},
    'H': {'3':'His','full':'Histidine'},
    'I': {'3':'Ile','full':'Isoleucine'},
    'L': {'3':'Leu','full':'Leucine'},
    'K': {'3':'Lys','full':'Lysine'},
    'M': {'3':'Met','full':'Methionine'},
    'F': {'3':'Phe','full':'Phenylalanine'},
    'P': {'3':'Pro','full':'Proline'},
    'S': {'3':'Ser','full':'Serine'},
    'T': {'3':'Thr','full':'Threonine'},
    'W': {'3':'Trp','full':'Tryptophan'},
    'Y': {'3':'Tyr','full':'Tyrosine'},
    'V': {'3':'Val','full':'Valine'},
    'U': {'3':'Sec','full':'Selenocysteine'},
    'O': {'3':'Pyl','full':'Pyrrolysine'},
    'B': {'3':'Asx','full':'Asparagine'},
    'Z': {'3':'Glx','full':'Glutamine'},
    'J': {'3':'Xle','full':'Leucine'},
    'X': {'3':'Xaa','full':'Unknown'}
}
############################################################
# DEFINITIONS FOR DIFFERENT SEQUENCE DATABASES
############################################################
UNIPROT_ALL = 0
UNIREF100 = 1
UNIREF90 = 2
UNIPROT_ALL_OLD = 3
UNIREF100_OLD = 4
UNIREF90_OLD = 5
BLAST_NEW_VERSION = 0
BLAST_OLD_VERSION = 1
############################################################
# PATHS FOR SEQUENCE DATABASES
############################################################
SEQ_DBS = {
    UNIPROT_ALL: {
        'fasta': '/clusterfs/ohana/external/UniProt/2013-12-17/uniprot_all.fasta', 
        'blastable_db': '/clusterfs/ohana/external/UniProt/2013-12-17/blastdbs/protein',
        'name': 'All UniProt (Dec. 2013) with new BLAST tools'
    },    
    UNIREF100: {
        'fasta': '/clusterfs/ohana/external/UniProt/2013-12-17/uniref100.fasta', 
        'blastable_db': '/clusterfs/ohana/external/UniProt/2013-12-17/blastdbs/protein_ur100',
        'name': 'UniRef100 (Dec. 2013) with new BLAST tools'
    },
    UNIREF90: {
        'fasta': '/clusterfs/ohana/external/UniProt/2013-12-17/uniref90.fasta', 
        'blastable_db': '/clusterfs/ohana/external/UniProt/2013-12-17/blastdbs/protein_ur90',
        'name': 'UniRef90 (Dec. 2013) with new BLAST tools'
    },
    UNIPROT_ALL_OLD: {
        'fasta': '/clusterfs/ohana/external/UniProt/2013-12-17/uniprot_all.fasta', 
        'blastable_db': '/clusterfs/ohana/external/UniProt/2013-12-17/blastdbs/prot',
        'name': 'All UniProt (Dec. 2013) with old BLAST tools'
    },    
    UNIREF100_OLD: {
        'fasta': '/clusterfs/ohana/external/UniProt/2013-12-17/uniref100.fasta', 
        'blastable_db': '/clusterfs/ohana/external/UniProt/2013-12-17/blastdbs/prot_ur100',
        'name': 'UniRef100 (Dec. 2013) with old BLAST tools'
    },
    UNIREF90_OLD: {
        'fasta': '/clusterfs/ohana/external/UniProt/2013-12-17/uniref90.fasta', 
        'blastable_db': '/clusterfs/ohana/external/UniProt/2013-12-17/blastdbs/prot_ur90',
        'name': 'UniRef90 (Dec. 2013) with old BLAST tools'
    }
}
###########################################################
# PDB LIBRARY - base path where pdb structures are located.
###########################################################
PDB_LIBRARY_PATH = '/clusterfs/ohana/external/pdb/'
