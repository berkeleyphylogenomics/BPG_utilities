from django.template import Context, RequestContext, loader
from django.shortcuts import render_to_response
from django.http import QueryDict, HttpResponse, HttpResponseRedirect
from django.forms.util import ErrorList
from django.utils.safestring import mark_safe

from pfacts003.utils.json_data_objects import annotatedTree

from pfacts003.phog.views import get_parameter_from_GET

def phyloscope(request, family_accession):
    t = loader.get_template('phyloscope/phyloscope.html')
    c = RequestContext(request, {'json_tree': annotatedTree(family_accession),
                                 'family_accession': family_accession,
                                 'hide_ortho_controls': True})
    return HttpResponse(t.render(c))

def phyloscope_get(request):
    t = loader.get_template('phyloscope/phyloscope.html')
    mode = get_parameter_from_GET('mode', request)
    subtreeLeftId = get_parameter_from_GET('subtreeLeftId', request)
    highlightLeftIds = get_parameter_from_GET('highlightLeftIds', request)
    colorLeftIds = get_parameter_from_GET('colorLeftIds', request)
    phog_accession =  get_parameter_from_GET('phog', request)
    family_accession =  get_parameter_from_GET('family', request)
    preferred_name = get_parameter_from_GET('preferredName', request)
    url_tree_method = get_parameter_from_GET('treeMethod', request)
    c = RequestContext(request, {
        'json_tree': annotatedTree(family_accession, url_tree_method),
        'family_accession': family_accession,
        'preferred_name': preferred_name,
        'url_tree_method': url_tree_method,
        'hide_ortho_controls': True})
    return HttpResponse(t.render(c))

def phyloscope_quick_start(request):
    t = loader.get_template('phyloscope/quick_start.html')
    c = RequestContext(request, {})
    return HttpResponse(t.render(c))

