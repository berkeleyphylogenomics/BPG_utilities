# Create your views here.

from pfacts003.phylofacts.models import UniProtTaxonomy
from pfacts003.phog.orthologs import OrthologTypes

from django.shortcuts import render_to_response

from pfacts003.phog.orthologs import getOrthologQuerySet
from pfacts003.phylofacts.models import PHOGRow

def index(request):
    
    return render_to_response('hpy/index.html', {
        'taxon': UniProtTaxonomy.objects.get(mnemonic='HELPY'),
        'js': {'jquery': True, 'expandable': True}
    })

def about(request):
    taxon = UniProtTaxonomy.objects.get(mnemonic='HELPY')
    return render_to_response('hpy/about.html', {'taxon': taxon})

def sequence(request):
    
    if 'sequence_id' in request.GET:

        uniprot, phogs, best_phogs, orthologs, phog_of_ortholog, error = \
            getOrthologQuerySet(
                request.GET['sequence_id'],
                request.GET.get('ortholog_type', OrthologTypes.PHOG_T_Medium),
                request.GET.get('threshold', 0.296874),
            )

        return render_to_response('hpy/sequence.html', {
            'phog_rows': [PHOGRow(phog) for phog in best_phogs],
        })

    return render_to_response('hpy/sequence.html', {})
