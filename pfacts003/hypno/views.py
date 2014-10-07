# Create your views here.
import os, urllib
from django.http import HttpResponse, Http404
from django.shortcuts import render_to_response
from django.template import RequestContext

def home(request):
    return render_to_response('hypno/HYPNO.htm', {}, context_instance=RequestContext(request))

def help(request):
    return render_to_response('hypno/HYPNO_help.htm', {}, context_instance=RequestContext(request))

def serve_hypno_source_code(request):
    response = HttpResponse(mimetype="application/zip")
    response['Content-Disposition'] = 'attachment; filename=HYPNO-master.zip'
    response.write(file('/clusterfs/ohana/software/hypno/data/HYPNO-master.zip', 'rb').read())
    return response

def serve_hypno_supplementary_data(request):
    response = HttpResponse(mimetype="application/x-compressed")
    response['Content-Disposition'] = 'attachment; filename=HYPNO_data.tgz'
    response.write(file('/clusterfs/ohana/software/hypno/data/HYPNO_data.tgz', 'rb').read())
    return response
