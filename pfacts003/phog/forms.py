import re

from django import forms
from django.forms.widgets import HiddenInput
from django.utils.safestring import mark_safe


from pfacts003.phylofacts.models import OrthologTypes, preset_thresholds, \
    get_ortholog_type_of_threshold
from pfacts003.utils.id_patterns import *
from pfacts003.utils.messages import unrecognized_identifier

OrthologChoices = [(OrthologTypes.SuperOrtholog, "Super orthology: 0.0"),
                    (OrthologTypes.PHOG_T_Tight, "Close: 0.09375"),
                    (OrthologTypes.PHOG_T_Medium, "Moderate: 0.296874"),
                    (OrthologTypes.PHOG_T_Loose, "Distant: 0.9375"),
                    (OrthologTypes.PHOG_T_Custom, "Custom: Enter value"),
                  ]

TAXON_CHOICES = [
      ('7165', 'Anopheles_gambiae'),
     ('63363', 'Aquifex_aeolicus'),
      ('3702', 'Arabidopsis_thaliana'),
    ('451804', 'Aspergillus_fumigatus'),
      ('1423', 'Bacillus_subtilis'),
       ('818', 'Bacteroides_thetaiotaomicron'),
       ('375', 'Bradyrhizobium_japonicum'),
      ('6239', 'Caenorhabditis_elegans'),
      ('5476', 'Candida_albicans'),
    ('315277', 'Chlamydia_trachomatis'),
    ('324602', 'Chloroflexus_aurantiacus'),
      ('1299', 'Deinococcus_radiodurans'),
    ('515635', 'Dictyoglomus_turgidum'),
      ('7227', 'Drosophila_melanogaster'),
     ('83333', 'Escherichia_coli'),
     ('76856', 'Fusobacterium_nucleatum'),
     ('35554', 'Geobacter_sulfurreducens'),
     ('33072', 'Gloeobacter_violaceus'),
    ('478009', 'Halobacterium_salinarum'),
      ('9606', 'Homo_sapiens'),
      ('6945', 'Ixodes_scapularis'),
    ('374847', 'Korarchaeum_cryptofilum'),
       ('173', 'Leptospira_interrogans'),
      ('2190', 'Methanocaldococcus_jannaschii'),
      ('2214', 'Methanosarcina_acetivorans'),
     ('10090', 'Mus_musculus'),
      ('1773', 'Mycobacterium_tuberculosis'),
    ('381754', 'Pseudomonas_aeruginosa'),
      ('4932', 'Saccharomyces_cerevisiae'),
      ('4896', 'Schizosaccharomyces_pombe'),
      ('1902', 'Streptomyces_coelicolor'),
      ('2287', 'Sulfolobus_solfataricus'),
      ('1148', 'Synechocystis_sp.'),
    ('289376', 'Thermodesulfovibrio_yellowstonii'),
      ('2336', 'Thermotoga_maritima'),
     ('69014', 'Thermococcus_kodakarensis'),
      ('5270', 'Ustilago_maydis'),
      ('4952', 'Yarrowia_lipolytica'),
  ]

class PairwiseGenomicOrthologyForm(forms.Form):
  taxon1 = forms.ChoiceField(choices=TAXON_CHOICES)
  taxon2 = forms.ChoiceField(choices=TAXON_CHOICES)
  def clean(self):
    taxon_name1 = self.cleaned_data['taxon1']
    taxon_name2 = self.cleaned_data['taxon2']
    if taxon_name1 == taxon_name2:
      raise forms.ValidationError('''Error: You must choose two distinct
        species for which to download the pairwise orthologs.''')
    return self.cleaned_data

class PHOGForm(forms.Form):
  phog_accession = forms.CharField(widget=forms.TextInput(attrs={'size':'40',
                                  'onfocus':"form.phog_accession.value='';"}),
                                  required=True, initial="PHOG accession")
  def clean_phog_accession(self):
    phog_accession = self.cleaned_data.get('phog_accession', '')
    if  uniprot_accession_re1.match(phog_accession) == None and \
        uniprot_accession_re2.match(phog_accession) == None and \
        uniprot_identifier_re1.match(phog_accession) == None and \
        uniprot_identifier_re2.match(phog_accession) == None and \
        uniprot_identifier_re3.match(phog_accession) == None and \
        gi_re.match(phog_accession) == None and \
        phog_re.match(phog_accession) == None and \
        scopid_re.match(phog_accession.lower()) == None and \
        bpgid_re.match(phog_accession.lower()) == None:
      raise forms.ValidationError(unrecognized_identifier)
    return phog_accession.upper()

class OrthologForm(forms.Form):
  ortholog_type = forms.ChoiceField(choices=OrthologChoices, required=False)
  ortholog_type.widget.attrs['onchange'] = mark_safe(
    "t=document.getElementById('id_threshold'); v=this.value; " 
    + "switch(v) {case '0': t.value='0.0'; break; case '1': t.value='0.09375'; " 
    + "break; case '2': t.value='0.296874'; break; case '3': t.value='0.9375'; " 
    + "break; case '4': t.value='Enter custom threshold'; break;}")

  sequence_id = forms.CharField(widget=forms.TextInput(attrs={'size':'40'}),
                                required=False)
  sequence_fasta = forms.CharField(widget=forms.Textarea(attrs={'cols':'50',
                                                                'rows':'6',
                                                                'style': 'font-size:10px'}),
                                    required=False)
  threshold = forms.FloatField(required=False, initial="0.0")
  threshold.widget.attrs['onclick'] = mark_safe("if (this.value=='Enter custom threshold') {this.value='';}")
  threshold.widget.attrs['onchange'] = mark_safe("document.getElementById('id_ortholog_type').selectedIndex=4;")

  def clean(self):
    if ('sequence_id' not in self.cleaned_data or \
        self.cleaned_data['sequence_id'] == u'') \
        and ('sequence_fasta' not in self.cleaned_data or \
            self.cleaned_data['sequence_fasta'] == u''):
      raise forms.ValidationError("Must enter either an identifier or a sequence in FASTA format")
    else:
      return self.cleaned_data
  def clean_ortholog_type(self):
    ortholog_type = self.cleaned_data.get('ortholog_type', 
                                          OrthologTypes.SuperOrtholog)
    if ortholog_type in preset_thresholds:
      threshold = preset_thresholds[ortholog_type]
    elif ortholog_type != OrthologTypes.PHOG_T_Custom:
      orthoog_type = OrthologTypes.SuperOrtholog
    return ortholog_type
  def clean_threshold(self):
    threshold = self.cleaned_data.get('threshold', 0.0)
    if not threshold:
      threshold = 0.0
      ortholog_type = OrthologTypes.SuperOrtholog
      return threshold
    if threshold < 0.0:
      raise forms.ValidationError("Threshold must be non-negative")
    ortholog_type = get_ortholog_type_of_threshold(threshold)
    return threshold
  def clean_sequence_id(self):
    sequence_id = self.cleaned_data.get('sequence_id', '').strip()
    if sequence_id == '':
      return sequence_id
    if  uniprot_accession_re1.match(sequence_id) == None and \
        uniprot_accession_re2.match(sequence_id) == None and \
        uniprot_identifier_re1.match(sequence_id) == None and \
        uniprot_identifier_re2.match(sequence_id) == None and \
        uniprot_identifier_re3.match(sequence_id) == None and \
        gi_re.match(sequence_id) == None and \
        phog_re.match(sequence_id) == None and \
        scopid_re.match(sequence_id.lower()) == None and \
        bpgid_re.match(sequence_id.lower()) == None:
      raise forms.ValidationError(unrecognized_identifier)
    return sequence_id.upper()
  def clean_sequence_fasta(self):
    sequence_fasta = self.cleaned_data.get('sequence_fasta', u'')
    if sequence_fasta == u'':
      return sequence_fasta
    if fasta_re.match(sequence_fasta) == None:
      raise forms.ValidationError("Invalid FASTA input")
    return sequence_fasta
