
from django import forms

from pfacts003.common.forms import fasta_form_class
from pfacts003.utils.format_patterns import format_protein_fasta

MIN_SEQUENCES = 2
MAX_SEQUENCES = 200

class IntrepidForm(fasta_form_class(name='intrepid_fasta',
    vargs={'min_seqs': MIN_SEQUENCES, 'max_seqs': MAX_SEQUENCES})):
    pass
