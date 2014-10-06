#!/usr/bin/env python

"""A quick and dirty script to validate Flower Power Results

    This script was created quickly and, although useful, will probably need
updating upon the next iteration of Quality Assurance (QA) of flower power
results. 

    The input to this file (DEFAULT_INPUT_FILE) describes seeds for each
subfamily, relative to the PFAM_ROOT. For example, the following three entries
are for PFAM accession PF00575.15. Three subfamilies are described (0, 1, and
10) and the seed sequence to that family is given in the file mentioned:

PF00575.15/subfam0/NUSA_THET8_141-208.fa
PF00575.15/subfam1/RNG_SHIFL_35-122.fa
PF00575.15/subfam10/RPOE_SULAC_78-150.fa

    For each of the above mentioned subfamilies, the final align-to-model (a2m) file is exampled. The filename for the a2m file is found by
get_a2m_filename(path). We review these a2m files checking:
    * For extremely long gaps,
    * For overall gappiness, 
    * For singletons, 
    * etc.

   This file can be improved by:
    - Using the logging function,
    - Using regular expressions where it may be more appropriate,
    - Fixing PEP-8 (http://www.python.org/dev/peps/pep-0008/) violations
    - Fixing redundant program logic
    - Better printing/error reporting

    The definitions of the QA script changed over time, allowing for confusing
logic left in this script. However, due to time contraints, we much let this go
until this is reviewed again.
"""

from optparse import OptionParser
import os
import re
import sys
import subprocess

from Bio import SeqIO


PFAM_ROOT = '/clusterfs/ohana/bpg/pfam'
DEFAULT_INPUT_FILE =\
    '/tmp/ordered_pfam_subfam_subdirs.txt'
FINAL_A2M = 'final.a2m'
GAP_LENGTH = 40 # Gap string is '-'*GAP_LENGTH
GAP_COUNT = 10


def get_a2m_filename(path):
    """Return full absoluate path to the align-to-model file"""

    (head, tail) = os.path.split(path)
    a2m_file = os.path.join(head, FINAL_A2M)

    return a2m_file


def skip_or_match(char):
    """Let uppercase letters and dashes be returned

    Hidden Markov Model profiles consist of upper case letters (matches), lower
case letters (insert states), dashes, and dots. This function will filter
(return) only the skip state (dash) and capital letters.
"""
    if char == '-':
        return char
    if char.isupper():
            return char

def get_sequence_length(filename):
    """Return the length of the first sequence in the filename given"""

    seed_file = open(filename, 'r')
    seed_sequence = SeqIO.parse(seed_file, "fasta")

    # Convert sequence to string
    first_sequence = '%s' % seed_sequence.next().seq

    # Get length of only Uppercase letters and gaps
    sequence_length = len(filter(skip_or_match, first_sequence))

    seed_file.close()

    return sequence_length


def get_original_length(hmm):
    """Return original sequence length from given hmm

    The sequence (hmm) passed can contain dots, dashes (skip states), upper
case letters (matches), and lower case letters (insert states). The original
length can be determined by looking at only upper case letters, and dashes.
This is what is returned from return_sequence_length.
"""

    test_length_string = filter(skip_or_match, hmm)
    num_gaps = test_length_string.count('-')

    return len(test_length_string), num_gaps


def iterate_sequences(filename, original_length, gap_string):
    """Iterate through sequences in file gathering statistics

    We are mostly interested in singletons (a file with a single sequence), and
with sequences that have gaps. We, will therefore, count both as we iterate
through the list.
"""

    count = 0
    gap_count = 0
    possible_gaps = 0
    long_gap_count = 0
    percent_gaps = 0.0
    my_fasta = open(filename, 'r')

    for record in SeqIO.parse(my_fasta, "fasta"):
        count = count + 1

        sequence_str = '%s' % record.seq

        # Check sequence length (only matches and gaps)
        sequence_length, line_gap_count = get_original_length(sequence_str)
        if not sequence_length == original_length:
            print "Error: sequence_length: %d should be %d. See sequence '%s' in file '%s'" %\
                (sequence_length, original_length, record.id, filename)

        possible_gaps = possible_gaps + sequence_length
        gap_count = gap_count + line_gap_count

        gap_site = sequence_str.find(gap_string)
        if gap_site > 0:
            long_gap_count = long_gap_count + 1
    my_fasta.close()

    percent_gap = (gap_count * 1.0/possible_gaps)*100
    return (count, long_gap_count, percent_gap)


def start_qa(input_file, gap_length, max_gaps, verbose=False):

    assert(isinstance(gap_length, int))

    print "Starting QA test of FlowerPower job with:"
    print "gap_length:", gap_length
    print "max_gaps:", max_gaps

    gap_string = '-'*gap_length

    pfam_open_list= open(input_file, 'r')
    for pfam in pfam_open_list:
        my_path = os.path.join(PFAM_ROOT, pfam.strip())
        a2m_file = get_a2m_filename(my_path)

        if not os.path.isfile(my_path):
            print "ERROR: fa file expected, but not found:", my_path

        if not os.path.isfile(a2m_file):
            print "ERROR: Final a2m file expected, but not found:", a2m_file
            continue

        sequence_length = get_sequence_length(my_path)
        (count, gap_count, gap_percentage) = iterate_sequences(a2m_file,
            sequence_length, gap_string)

        if count < 2:
            print "SINGLETON: ", a2m_file
        if gap_count > max_gaps:
            print "GAPPY: ", a2m_file
        if verbose == True:
            print "%d\t%d\t%5.2f%%\t%s" % (count, gap_count, gap_percentage,
                a2m_file)

    print "End QA Test..."


def main():
    """Options are passed and start_qa is called"""

    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-i', '--input', dest='input_file',
        help="The input file indicating pfam data to check",
        default=DEFAULT_INPUT_FILE)
    parser.add_option('-v', '--verbose', dest='verbose',
        action="store_true",
        help="For more verbose output",
        default=False)
    parser.add_option('-l', '--gap-length', dest='gap_length',
        help="The length of a string of dashes('-') that is considered gappy",
        default=GAP_LENGTH)
    parser.add_option('-c', '--gap-count', dest='gap_count',
        help="The number of 'gaps' before a subfamily becomes 'gappy'",
        default=GAP_COUNT)

    (options, args) = parser.parse_args()

    if options.verbose:
        verbose = True
    else:
        verbose = False

    try:
        gap_length = int(options.gap_length)
    except ValueError:
        print "Error: Gap length should be an integer. However, found '%s'" % options.gap_length
        sys.exit(1)

    try:
        gap_count = int(options.gap_count)
    except ValueError:
        print "Error: Gap Count should be an integer. However, found '%s'" % options.gap_count
        sys.exit(1)

    ret = start_qa(options.input_file, gap_length, gap_count, verbose)


if __name__ == '__main__':
    main()
