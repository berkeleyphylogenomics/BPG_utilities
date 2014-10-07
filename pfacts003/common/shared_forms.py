from django import forms
from django.forms import ValidationError
from django.utils.safestring import SafeString
import os, re

fasta_seq_pat   = re.compile(r'^([AC-IK-NP-TVWY\s]+)$')
fasta_junk_pat  = re.compile(r'[\s.-]')
pdbid_chain_pat = re.compile(r'^.{4}:[a-z]$')

"""
text-area field that handles sequence data
"""
class SequenceTextArea(forms.CharField):

    def __init__(self):
        forms.CharField.__init__(self, min_length=1,
                                 widget=forms.widgets.Textarea({'cols': '100', 'rows': '10'}))

    # overrides the clean function to make sure the
    # sequence is in the proper format
    def clean(self, value):
        
        newval = value.strip()

        # Split the header and sequence (define a default header if necessary)
        if newval and newval.startswith('>') and newval.count('\n'):
            header, sequence = newval.replace('\r','\n').split('\n', 1)
        else:
            header, sequence = ('>query', newval)

        # Remove junk characters from the sequence and capitalize
        sequence = fasta_junk_pat.sub('', sequence).upper()

        # Raise exception for an empty sequence
        if not sequence:
            raise ValidationError('A protein sequence is required')

        # Raise exception for an invalid sequence
        elif not fasta_seq_pat.match(sequence):
            raise ValidationError(
                'Putative protein sequence has illegal residue codes')
        
        # Raise exception for a suspected nucleotide sequence
        elif float( sequence.count('A') + \
                    sequence.count('C') + \
                    sequence.count('G') + \
                    sequence.count('T') + \
                    sequence.count('U') )/len(sequence) > 0.9:
            raise ValidationError('Putative protein sequence appears to be nucleotide')

        return '%s\n%s\n' % (header,sequence)


"""
text-area field that handles PDB identifier data
"""
class PDBCharField(forms.CharField):

    def __init__(self):
        forms.CharField.__init__(self, min_length=6, max_length=6)

    # overrides the clean function to make sure the
    # sequence is in the proper format
    def clean(self, value):
        newval = value.strip().lower()

        # empty field is an error
        if not pdbid_chain_pat.match(newval):
            if not newval:
                raise ValidationError('A valid pdb entry is required')
            raise ValidationError('PDB entry must by of the form: PDBID:CHAIN')

        newval = newval[0:-1] + newval[-1].upper()

        pdbid   = newval[0:4] 
        chainid = newval[5]

        # check that the pdb file is there
        pdbfile = '/clusterfs/ohana/external/pdb_entries/%s/pdb%s.ent' % (pdbid[1:3], pdbid)

        if not os.path.exists(pdbfile):
            raise ValidationError('We\'re sorry, we don\'t appear to have a structure for pdb identifier %s' % pdbid)

        return newval


"""
simple input form that has a (protein) sequence in FASTA and an
(optional) email address for returning results
"""
class SequenceInputForm(forms.Form):
    sequence  = SequenceTextArea()
    useremail = forms.EmailField(required=False)
    workdir   = forms.CharField(required=False)

    def __init__(self, *arglist):
        forms.Form.__init__(self,*arglist)

        self.seq_field_keys = ('sequence','useremail','workdir')

    def as_hidden(self):
        mystr = ''

        mystr += self['sequence'].as_hidden() + '\n'
        mystr += self['useremail'].as_hidden() + '\n'
        mystr += self['workdir'].as_hidden() + '\n'

        return SafeString(mystr)

    def as_table(self):

        mytable = ''

        # print errors
        if self._errors:

            mytable += '<tr><td colspan="2" style="color: red;">\n'

            for errkey in self._errors:

                mytable += '%s:%s\n' % (self[errkey].label, self._errors[errkey])

            mytable += '</td></tr>'

        mytable += \
        """
        <tr>
            <td colspan="2">
                %s
                <b>Paste protein sequence in
                <a href="http://en.wikipedia.org/wiki/Fasta_format">FASTA format</a></b>
            </td>
        </tr>

        <tr>
            <td colspan="2">
                %s
            </td>
        </tr>
        <tr>
            <td>
                %s&nbsp;<label for="id_useremail">email results to you (optional)</label>
            </td>
            <td>
                <input style="width: 6em;" type="submit" value="submit"\>
                <input style="width: 6em;" type="reset" value="reset"\>
            </td>
        </tr>
        """ % (self['workdir'].as_hidden(), self['sequence'], self['useremail'])

        return SafeString(mytable)


""" 
simple input form that has a PDB identifier and an (optional)
email address
"""
class PDBInputForm(SequenceInputForm):
    sequence  = PDBCharField()
    useremail = forms.EmailField(required=False)
    workdir   = forms.CharField(required=False)

    def as_table(self):

        mytable = ''

        # get any errors
        if self._errors:
            if 'sequence' in self._errors:
                mytable += '<tr><td colspan=3 style="color: red;">ERROR: %s</td></tr>' % self._errors['sequence']

            if 'usermail' in self._errors:
                mytable += '<tr><td colspan=3 style="color: red;">ERROR: %s</td></tr>' % self._errors['useremail']

        mytable += \
        """
        <tr>
            <td>
                %s
                <b>Input PDB sequence as PDBID:CHAIN:</b>
                <div style="font-size: 10pt;">
                For example, <q>12as:A</q> will look up <br>
                chain A for PDB entry 12as
                </div>
            </td>
            <td style="vertical-align: top;">
                %s
            </td>
            <td>
                <input style="width: 6em;" type="submit" value="submit"\>
                <br>
                <input style="width: 6em;" type="reset" value="reset"\>
            </td>
        </tr>
        <tr>
            <td></td>
        </tr>
        <tr>
            <td colspan=3>
                %s email results to you (optional)
            </td>
        </tr>
        """ % (self['workdir'].as_hidden(), self['sequence'], self['useremail'])

        return SafeString(mytable)

