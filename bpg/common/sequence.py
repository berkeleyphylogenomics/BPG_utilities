#!/usr/bin/python
import re, hashlib, base64
from pfacts003.phylofacts.models import Family

# ------------------------------------------------------------------------------
def is_dna_or_rna(seq):
    """Verifies wether a seq is DNA or RNA (composed only of [acgtuxn] chars).
       Inputs:
       seq = query sequence.
       Outputs:
       is_dna_rna_f = boolean indicating whether sequence is dna or rna.
    """
    is_dna_rna_f = False
    if re.compile('[^a-zA-Z]+').match(seq) :
        print seq
        raise ValueError, 'Sequence argument <seq> has invalid characters.'
    if len(seq) < 1 :
        raise ValueError, 'Invalid empty sequence!'
    if re.compile('^[acgtuxnACGTUXN]+$').match(seq):
        is_dna_rna_f = True
    else:
        is_dna_rna_f = False
    return is_dna_rna_f

# ------------------------------------------------------------------------------
def is_mostly_dna(seq, ratio):
    """Verifies wether a seq is mostly DNA or RNA (composed mostly of [acgtu] chars).
       Inputs:
       seq = query sequence.
       ratio = min ratio between dna_content and seq_length that will cause seq
               to be flagged as dna.
       Outputs:
       is_mostly_dna_f = boolean indicating whether sequence is mostly dna or rna.
    """
    is_mostly_dna_f = False
    if ratio < 0 or ratio > 1 :
        print 'ratio = %d' % ratio
        raise ValueError, '<ratio> has to be between 0 and 1.'
    if re.compile('[^a-zA-Z]+').match(seq) :
        print seq
        raise ValueError, 'Sequence argument <seq> has invalid characters.'
    if len(seq) < 1 :
        raise ValueError, 'Invalid empty sequence!'
    dna_letters = 0
    for letter in seq :
        if re.compile('[acgtuACGTU]').match(letter) :
            dna_letters = dna_letters + 1
    final_ratio = float(dna_letters)/float(len(seq))
    if final_ratio > ratio :
        is_mostly_dna_f = True
    return is_mostly_dna_f

# -----------------------------------------------------------------

# -----------------------------------------------------------------
def get_family_of_sequence_by_hash(seq):
# This function hashes the input sequence seq with a sha1 hash and encodes it in
# base 64.  The return is a list of distinct family accessions that the sequence is in.

    hashobj = hashlib.sha1()
    hashobj.update(seq)
    seqhash = base64.b64encode(hashobj.digest()).strip('=')
    return [x.get_accession() for x in Family.objects.filter(trees__treenodes__sequence_header__sequence__seguid=seqhash).distinct()]
    
