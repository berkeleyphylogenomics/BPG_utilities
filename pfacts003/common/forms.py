"""Reusable Form Fields with Validation"""

from copy import deepcopy

from django import forms

from pfacts003.utils.format_patterns import FastaError, format_protein_fasta

def fasta_form_class(name, vargs=None, upload=True):
    """Generic FASTA Form with validation from optional validation arguments

    Validation Arguments (vargs) is passed to FastaField. The following
    parameters may be used:

    * to_upper (default=False) Converts sequences to upper case letters
    * wrap (default=False) Reset the line wrapping for sequences
    * alignment (default=False) Allow '-' character and enforce equal lengths
    * allow (default=None) String of additional characters allowed in the sequence
    * extra_chars (default=False) Allow BOUXZ
    * stops (default=False) Allow '*' character
    * ignore_case (default=False) Check sequences insensitive to case
    * min_seqs (default=1) Minimum number of sequences allowed
    * max_seqs (default=1) Maximum number of sequences allowed
    """

    def fasta_clean(self):
        cl = type(self)

        ft, ff, fe, va = getattr(cl, '_fasta_handler')
        
        cdata = self.cleaned_data
        ftext = cdata.get(ft, False)
        ffile = cdata.get(ff, False)

        if not bool(ftext) ^ bool(ffile):
            raise forms.ValidationError('FASTA data required; submit either as text or file input, not both.')

        try:
            if bool(ftext):
                cdata[name] = format_protein_fasta(ftext, **va)
                del cdata[ft]
            else:
                tmp = ffile.read()
                if getattr(ff, 'delete', False):
                    ffile.delete(save=False)
                del cdata[ff]
                cdata[name] = format_protein_fasta(tmp, **va)
        except FastaError, err:
            raise forms.ValidationError("Invalid FASTA input (%s; see %s)" %\
                (err.message, err.line))

        return cdata

    fasta_text = name + '_fasta_text'
    fasta_file = name + '_fasta_file'
    if vargs is None:
        vargs = {}

    attrs = {
        fasta_text: forms.CharField(required=False, label='Paste FASTA',
            widget=forms.Textarea(attrs={'cols':'70', 'rows':'6'})),
        fasta_file: forms.FileField(required=False, label='Upload FASTA',
            widget=upload and forms.FileInput or forms.HiddenInput),
        '_fasta_handler': (
            # the name of the text field
            name + '_fasta_text',
            # the name of the file field
            name + '_fasta_file',
            # the name of the entry in cleaned_data
            name,
            # pass blank dict if vargs is None, otherwise pass vargs
            vargs is None and {} or vargs
        ),
        'clean': fasta_clean,
        '_fasta_clean': fasta_clean,
    }

    return type('FastaForm', (forms.Form,), attrs)

SingleProteinForm = fasta_form_class('protein', vargs={
        'min_seqs': 1,
        'max_seqs': 1,
    })

AlignToModelForm = fasta_form_class('a2m', vargs={
        'alignment': True,
        'allow': '.',
        'min_seqs': 2,
        'max_seqs': 50,
    })

class ContactForm(forms.Form):
    sender = forms.EmailField(required=False)
    cc_myself = forms.BooleanField(required=False)
    subject = forms.CharField(max_length=100, required=False)
    message = forms.CharField(
        widget=forms.Textarea(attrs={'cols':'63', 'rows':'6'},),
        required=True)

    def clean(self):
        cleaned_data = self.cleaned_data
        sender = cleaned_data.get("sender")
        cc_myself = cleaned_data.get("cc_myself")

        if cc_myself and not sender:
            raise forms.ValidationError(\
                "We can't cc you without an email address")

        subject = cleaned_data.get("subject")
        if not cleaned_data['subject']:
            cleaned_data['subject'] = "BPG 'Contact Us' Email"

        return cleaned_data
