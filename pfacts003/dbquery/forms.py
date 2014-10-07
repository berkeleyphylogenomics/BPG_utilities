from django import forms
from urllib import urlencode
from pfacts003.phylofacts.models import TreeNode, UniProtTaxonomy, GO_Evidence, UniProtGO, UniProtEC
from django.db.models import Q
from dbquery_results import GetGenesInTaxonWithOrthologsWithGOAnnotation, GetGenesInTaxonWithOrthologsWithEC

_evidence_codes = (
    # Experimental
    ('EXP', 'EXP: Inferred from Experiment'),
    ('IDA', 'IDA: Inferred from Direct Assay'),
    ('IPI', 'IPI: Inferred from Physical Interaction'),
    ('IMP', 'IMP: Inferred from Mutant Phenotype'),
    ('IGI', 'IGI: Inferred from Genetic Interaction'),
    ('IEP', 'IEP: Inferred from Expression Pattern'),
    # Computational Analysis
    ('ISS', 'ISS: Inferred from Sequence or Structural Similarity'),
    ('ISO', 'ISO: Inferred from Sequence Orthology'),
    ('ISA', 'ISA: Inferred from Sequence Alignment'),
    ('ISM', 'ISM: Inferred from Sequence Model'),
    ('IGC', 'IGC: Inferred from Genomic Context'),
    ('RCA', 'RCA: inferred from Reviewed Computational Analysis'),
    # Author Statement
    ('TAS', 'TAS: Traceable Author Statement'),
    ('NAS', 'NAS: Non-traceable Author Statement'),
    # Curator Statement
    ('IC', 'IC: Inferred by Curator'),
    ('ND', 'ND: No biological Data available'),
    # Automatically-assigne
    ('IEA', 'IEA: Inferred from Electronic Annotation'),
    # Obsolete
    ('NR', 'NR: Not Recorded'),
)

_evidence_choices = (
    ('Experimental Evidence Codes', _evidence_codes[0:6]),
    ('Computational Analysis Evidence Codes', _evidence_codes[6:12]),
    ('Author Statement Evidence Codes', _evidence_codes[12:14]),
    ('Curator Statement Evidence Codes', _evidence_codes[14:16]),
    ('Obsolete Evidence Codes', _evidence_codes[16:]),
    ('Other Options', (('ANY', 'Any Evidence'),)),
)

class _TaxonForm(forms.Form):
    taxon = forms.CharField()
    def as_url(self):
        return urlencode(self.cleaned_data)
    def _get_taxon(self):
        
        cd = self.cleaned_data
        tval = cd.get('taxon', False)    

        if not tval:
            del cd['taxon']
            return UniProtTaxonomy.objects.none()

        uptquery = Q(common_name__iexact=tval) | \
            Q(scientific_name__iexact=tval) | \
            Q(mnemonic__iexact=tval)
        try:
            uptquery |= Q(id=int(tval))
        except ValueError:
            pass
        uptquery = UniProtTaxonomy.objects.filter(uptquery)
        if not bool(uptquery):
            return UniProtTaxonomy.objects.none()

        #this really should be done in the db -- ideas?
        parents = []
        for item in uptquery.order_by('left_id', '-right_id').all():
            if not parents or parents[-1].right_id < item.right_id:
                parents.append(item)

        uptquery = UniProtTaxonomy.objects.filter(reduce(
            lambda x,y: x|y,
            map(
                lambda p: Q(
                    left_id__gte=p.left_id,
                    right_id__lte=p.right_id
                ),
                parents
            )
        ))

        if not bool(uptquery):
            self._errors['taxon'] = \
                self.error_class(["Could not find taxon '%s'." % cd['taxon']])
            del cd['taxon']
        
        return uptquery

class GoForm(_TaxonForm):
    term_type = forms.ChoiceField(label='GO Term Type', choices=(
        ('molecular_function', 'Molecular Function'),
        ('biological_process', 'Biological Process'),
        ('cellular_component', 'Cellular Localization'),
    ))
    name = forms.CharField(label='GO Term Name')
    evidence_code = forms.ChoiceField(choices=_evidence_choices,
        label='GO Evidence Code', initial='ANY')

    def clean(self):
        cd = self.cleaned_data
        for key in ('taxon','name','evidence_code','term_type'):
            if cd.has_key(key):
                cd[key] = cd[key].strip()

        # Find Taxon
        uptquery = self._get_taxon()

        # Find UniProtGO
        if cd.get('evidence_code', False) and cd.get('name', False) \
            and cd.get('term_type', False):
            evcode = cd['evidence_code']
            upgquery = UniProtGO.objects.filter(
                go_term__term_type=cd['term_type'],
                go_term__name__iexact=cd['name'],
            )
            if evcode == 'EXP':
                upgquery = upgquery.filter(go_evidence__id__lte=6)
            elif evcode != 'ANY':
                upgquery = upgquery.filter(go_evidence__evidence=evcode)
            if not bool(upgquery):
                self._errors['name'] = \
                    self.error_class(["Could not find term '%s' with evidence code '%s'." % (cd['name'], evcode)])
                del cd['name']
                del cd['evidence_code']
        if bool(self._errors):
            return cd

        self.results = GetGenesInTaxonWithOrthologsWithGOAnnotation(
            uptquery, upgquery,
        )
        if not bool(self.results[0]):
            raise forms.ValidationError('No sequences matched your query.')
        
        return cd

class ECForm(_TaxonForm):
    ec = forms.RegexField(r'^\d+(\.(\d+|-)){0,3}$', label='EC Number')
    brenda = forms.ChoiceField(label='BRENDA', choices=(
        ('False','Optional'),
        ('True','Required'),
    ))
    def clean(self):
        cd = self.cleaned_data
        for key in ('taxon','ec'):
            if cd.has_key(key):
                cd[key] = cd[key].strip()

        # Find Taxon
        uptquery = self._get_taxon()
        
        if cd.get('brenda', False) and cd.get('ec', False):
            upequery = UniProtEC.objects
            if cd['brenda'] == 'True':
                upequery = upequery.filter(is_in_brenda_f=True)
            upequery = upequery.filter(**dict(
                ('ec__%s_number' % f, n) for f,n in zip(
                    ('class', 'subclass', 'subsubclass', 'enzyme'),
                    cd['ec'].rstrip('.-').split('.'),
                )
            ))
            if not bool(upequery):
                self._errors['ec'] = \
                    self.error_class(["Could not find ec '%s'." % cd['ec']])
                del cd['ec']

        if bool(self._errors):
            return cd
        
        self.results = GetGenesInTaxonWithOrthologsWithEC(uptquery, upequery)
        if not bool(self.results[0]):
            raise forms.ValidationError('No sequences matched your query.')
        
        return cd
