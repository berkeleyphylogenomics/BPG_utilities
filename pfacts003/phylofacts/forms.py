from django import forms
import re
from pfacts003.utils.id_patterns import is_uniprot_identifier_format, is_uniprot_accession_format
from pfacts003.utils import messages
from pfacts003.common.forms import fasta_form_class, SingleProteinForm
from pfacts003.utils.extract_sequence_from_fasta import extract_sequence_from_fasta

class SequenceSearchForm(fasta_form_class(
        name='sequence_search', upload=False,
        vargs={
            'min_seqs': 1, 'bad_header_chars': '\(\)', 'max_seqs': 1,
            'to_upper': True, 'max_length': 2000,
        })):
    pass

class NewSequenceSearchForm(forms.Form):
    input = forms.CharField(widget = forms.Textarea);

    def clean_input(self):
        input = self.cleaned_data['input']
        defline, sequence, errors = extract_sequence_from_fasta(input)
        if sequence:
            if set([char.upper() for char in sequence]).issubset(set(['A', 'C', 'T', 'G'])):
                raise forms.ValidationError('''Error: You appear to have submitted a nucleotide 
                    sequence. Only amino acid sequences are accepted.''')
        else:
            if not errors:
                raise forms.ValidationError("You did not submit a sequence")
            error_list = []
            for error in errors:
                error_list.append("""position %d,('%s')""" % (error[1] + 1, error[0]))  
            error_list = ', '.join(error_list)
            if ">" in error_list:
                raise forms.ValidationError('''Error: illegal characters were found at: %s
                    . \nYou may have submitted more than one sequence.''' % error_list)
                raise forms.ValidationError("Error: illegal characters were found at: %s" % error_list)
            raise forms.ValidationError("Error: illegal characters were found at: %s" % error_list)
        return input


class FamilyAccessionForm(forms.Form):
    family_accession = forms.CharField()
    
    def clean_family_accession(self):
        family_accession = self.cleaned_data['family_accession'].strip()
        pattern = re.compile(r"^bpg\d{7}$", re.IGNORECASE) 
        if pattern.match(family_accession) is None:
            raise forms.ValidationError("Error: %s" % messages.family_accession_format_message)
        return family_accession 

class UniProtForm(forms.Form):
    acc_or_ident = forms.CharField()

    def clean_acc_or_ident(self):
      acc_or_ident = self.cleaned_data['acc_or_ident'].strip()
      print "in UniProtForm..."
      print "acc_or_ident = %s" % acc_or_ident
      if not (is_uniprot_identifier_format(acc_or_ident) or is_uniprot_accession_format(acc_or_ident)):
          raise forms.ValidationError("Error: %s" % messages.invalid_acc_or_ident_message)
      return acc_or_ident
 
