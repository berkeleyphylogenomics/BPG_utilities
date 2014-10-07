from django import forms

from pfacts003.common.forms import fasta_form_class
from pfacts003.utils.format_patterns import format_protein_fasta

class SatchmoForm(fasta_form_class(
        name='satchmo_fasta', upload=False,
        vargs={
            'min_seqs': 4, 'bad_header_chars': '\(\)', 'max_seqs': 500,
            'allow': 'XB', 'to_upper': True, 'max_length': 2000,
        })):
    pass

class AdvancedForm(forms.Form):
    # the first half of the field name corresponds to a 'step' from queued.py;
    # the second half corresponds to the attribute put into func_info.fields
    # for that step
    satchmo__minaff = forms.FloatField(initial='-3.00',
        label='Minimum Average Affinity')
    quicktree__treebuilder = forms.ChoiceField(choices=(
        ('FastTree', 'FastTree'),
        ('quicktree', 'QuickTree'),
    ), initial='quicktree', label='Tree Building Method')
