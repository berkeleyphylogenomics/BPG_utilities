#!/usr/bin/env python

"""Verify families on filesystem

This module is used as a library (import), or a silent cron script only. For
cron script use assistance, execute this module with --help.
    
For interactive use, see the manage_ohana.py master script.
""" 

from datetime import datetime
from optparse import OptionParser
import os
import re
import sys

from Bio import SeqIO

try:
    from bpg.manage_ohana.shared import logging, press_enter
except ImportError:
    sys.stderr.write("Activate environment first.")
    sys.exit(1)
    

# HERE - change back when put in production
FAMILY_ROOT = '/clusterfs/ohana/bpg/pfacts'
#FAMILY_ROOT = '/home/glenjarvis/experiments/sample_data/pfacts'
BOOK_LIST_FILENAME = '/home/glenjarvis/experiments/sample_data/pfacts/family_list_production.txt'


# If the following subdirectories or files don't exist within a file, then
# an error is reported.
required_files = [
    "bpg%06d.a2m", # a2m == Align to Model (Family MSA)
    "bpg%06d.fa", # fa == FASTA (consensus sequence)
    "bpg%06d.pfam", # protein family info, result from running hmm_pfam on fa
    "bpg%06d.phobius", # Transmembrane Domain predictions
    "bpg%06d.pdb.selected.dist", # Result of scoring family HMM against PDB
    "bpg%06d.astats", # astats == Alignment Statistics abou the family MSA
    "bpg%06d.mod", # General HMMs (GMMs) for family (not same as subfamily HMMS)
    "bpg%06d.nj", # nj == Neighbor Joining; NJ Tree for family
    "bpg%06d_sw2.mlib", # sw == Smith-Waterman; Calibration of the HMM
    "shmms/bpg%06d.subfam", # sw == Smith-Waterman; Calibration of the HMM
    "shmms/bpg%06d.tree.gz", # Shmm trees
    "shmms/bpg%06d.sfreps", # Subfamily Consensus Sequences
]

warning_files = [
    "bpg%06d.seed", # Seed sequence used to build this family
    "bpg%06d.pfam_domain_map.png", # Graphic generated by bliss program
    "bpg%06d_original.subfam", # Original subfam file
    "shmms/bpg%06d.ndist", 
    "shmms/bpg%06d.gtree.gz", 
    "shmms/bpg%06d.xml.gz", 
]

# We are no longer interested in any of these files
ignore_files = [
    'shmms/bpg%d.cost',
    'shmms/bpg%d.fit',
    'bpg%06d_nr90.a2m',
    'bpg%06d_nr100.a2m',
    'bpg%06d_pfam_colorfile',
    'bpg%06d_bpgprotseq.a2m',
    'bpg%06d_pfam_imagemap',
]

def get_family_path(accession):
    """Get path to family, given an accession

    If the accession does not match the expected format (bpg######), then None
    is returned. Otherwise, a full path to family is generated from FAMILY_ROOT
    and the accession.
    """

    if not re.match('^bpg\d{6}$', accession):
        return None

    # Generate path/Path will be {BASE}/bpg{3#s}/bpg{6#s}/user
    return os.path.normpath("%s/bpg%s/%s" % (FAMILY_ROOT,
        accession[3:6], accession))

def get_a2m_filename(accession):

    """Return filename for align to model file"""

    return os.path.join(get_family_path(accession), "%s.a2m" % accession)
    

def check_a2m(accession):

    """Check Align-To-Model (a2m)

    Checks all alignments in Align-To-Model (a2m) file, to ensure they are at
    least syntatically correct enough to be read by BioPython.
    """

    a2m_filename = get_a2m_filename(accession)

    # Don't report an error if a2m file is missing. This is already done
    # in check_one_family()
    if os.path.exists(a2m_filename):
        a2m = open(a2m_filename, 'r')
        for record in SeqIO.parse(a2m, "fasta"):
            # Nothing is done with the record. This is just to ensure it can
            # be parsed and read. If this fails, we'll know what Exception
            # to put here and log
            dummy = record.id
        a2m.close()


def check_one_family(accession, verbosity=False):
    """Check an individual family for required files

    A list of file paths to check is generated from the FAMILY_ROOT constant,
    the file_list given, and manipulation of the accession. After the list is
    generated, each path is checked for existance. If the path does not exist,
    an error is logged and possibly written to standard out.

    The following assumptions are expected:
        * accession is a string in the format of bpg###### (e.g., bpg000643)
	* file_list is a list of strings and each string has a six-digit
	  decimal print formatting specifier (e.g., "bpg%06d.seed")
    """

    if verbosity:
        sys.stdout.write("Checking family %s\n" % accession)

    this_family_dir = get_family_path(accession)

    # Build list of required files for family
    this_family_files =\
        [os.path.join(this_family_dir, file % int(accession[3:]))\
                                                    for file in required_files]

    for  file in this_family_files: 
        if not os.path.exists(file):
            if verbosity:
                sys.stderr.write(\
                             "\tFamily is corrupt. Missing file '%s'\n" % file)
            logging.error(\
             '%s: Family is corrupt. Missing file %s' % (datetime.now(), file))

    # Check the Alignment file in the family
    check_a2m(accession)

def check_families(input_file, verbosity=True):
    """Read list of families and iterate call check_one_family on each"""

    logging.info("Beginning families check.")
    num_corrupt = 0
    family_list_file = open(input_file, 'r')
    for family in family_list_file:
        if not check_one_family(family.strip(), verbosity):
            num_corrupt = num_corrupt + 1

    logging.info(\
        "Finished checking families. %d families corrupt" % num_corrupt)
    if verbosity:
        sys.stdout.write("\n\t%d Families were corrupted.\n\n" % num_corrupt)
        press_enter()


def main():
    """Options are passed and get_credentials called"""
    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-s', '--scriptable', dest='scriptable',
        action="store_true",
        help="For cron job scripting - no user interaction is needed",
        default=False)
    parser.add_option('-i', '--input', dest='input_file',
        help="The input file of families to check")
    parser.add_option('-u', '--username', dest='username',
        help="The username you wish to get credentials for")

    (options, args) = parser.parse_args()

    if options.username is not None:
        password = get_credentials(options.username)


    if options.input_file is None:
        input_file = BOOK_LIST_FILENAME
    else:
        input_file = options.input_file

    if not options.scriptable:
        print __doc__
    else:
        check_families(input_file, verbosity=False)


if __name__ == '__main__':
    main()
