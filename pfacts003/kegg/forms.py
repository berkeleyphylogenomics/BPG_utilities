from django import forms
from pfacts003.phog.forms import OrthologChoices
from pfacts003.utils.id_patterns import taxon_id_re
from django.utils.safestring import mark_safe
from pfacts003.phylofacts.models import OrthologTypes, preset_thresholds, \
    get_ortholog_type_of_threshold
from pfacts003.utils.messages import unrecognized_identifier

class KeggMapForm(forms.Form):
  require_brenda = forms.BooleanField(initial=True,required=False)
  ortholog_type = forms.ChoiceField(choices=OrthologChoices, required=False)
  ortholog_type.widget.attrs['onchange'] = mark_safe(
    "t=document.getElementById('id_threshold'); v=this.value; " 
    + "switch(v) {case '0': t.value='0.0'; break; case '1': t.value='0.09375'; " 
    + "break; case '2': t.value='0.296874'; break; case '3': t.value='0.9375'; " 
    + "break; case '4': t.value='Enter custom threshold'; break;}")

  taxon_id = forms.CharField(widget=forms.TextInput(attrs={'size':'40'}),
                                required=True)
  threshold = forms.FloatField(required=False, initial="0.0")
  threshold.widget.attrs['onclick'] = mark_safe("if (this.value=='Enter custom threshold') {this.value='';}")
  threshold.widget.attrs['onchange'] = mark_safe("document.getElementById('id_ortholog_type').selectedIndex=4;")

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
  def clean_taxon_id(self):
    taxon_id = self.cleaned_data.get('taxon_id', '').strip()
    if taxon_id == '':
      return taxon_id
    if taxon_id_re.match(taxon_id) == None:
      raise forms.ValidationError(unrecognized_identifier)
    return taxon_id
