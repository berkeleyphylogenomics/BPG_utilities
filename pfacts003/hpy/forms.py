from django import forms
from pfacts003.phog.orthologs import OrthologTypes

class SequenceForm(forms.Form):
    seqid = forms.RegexField(regex=r'^\w+(_\w+)?$', label='UniProt ID')
    threshold = forms.FloatField(widget=forms.HiddenInput,
        initial=0.296874)
    ortholog_type = forms.IntegerField(widget=forms.HiddenInput,
        initial=OrthologTypes.PHOG_T_Medium)
