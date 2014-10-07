#!/usr/bin/python
import re
import unittest

def extract_sequence_from_fasta(input, allowed_chars = 'ACDEFGHIKLMNPQRSTVWYX'):
    """
    Function takes input that may contain an amino-acid sequence, runs validation
    and tries to extract a defline, a useable sequence, and any errors 
    associated with any problematic inputs.  

    See unittests for this function for details on expected behavior with 
    different inputs.
    """
    lines = input.split("\n")
    first_line, remainder = lines[0], lines[1:]
    if first_line.strip() and first_line.strip()[0] == '>':
        defline = first_line.strip()
        sequence = ''.join(remainder)
    else:
        defline = ''
        sequence = ''.join([first_line, ''.join(remainder)])

    p = re.compile(r"[\s\d]*")
    sequence = p.sub('',sequence)

    # Allow trailing "*"
    if sequence and sequence[-1] == "*":
        sequence = sequence[:-1]

    allowed_chars = set([char for char in allowed_chars.upper()])
    sequence_length = len(sequence)
    errors = []
    for i in xrange(sequence_length):
        if sequence[i].upper() not in allowed_chars:
            errors.append((sequence[i], i))

    if errors:
        sequence = ''

    return (defline, sequence, errors)


