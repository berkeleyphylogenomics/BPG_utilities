# Create your views here.
from django.http import Http404, HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
from pfacts003.phylo4j.graph_queries import get_interactions, get_all_relationships
from pfacts003.phylo4j.consts import *
from textwrap import wrap
from StringIO import StringIO
from Bio import Phylo
from py2neo import neo4j

def index(request):
    return render_to_response('phylo4j/index.html', {}, context_instance = RequestContext(request)) 

def ppi(request):
    getstr = request.GET['query'].strip()
    # get the intact interactions
    interactions = get_interactions(getstr)
    return render_to_response('phylo4j/ppi.html', {
            'interactions': interactions,
            'query':  getstr 
        }, context_instance = RequestContext(request))

def explorer(request):
    return render_to_response('phylo4j/explorer.html', {
            'node_id': id,
            'all_relationships': get_all_relationships()
        }, context_instance = RequestContext(request))
