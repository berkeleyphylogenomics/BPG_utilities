import os

from django.shortcuts import render_to_response

from pfacts003.utils.files import serve_file

def index_pqq(request):
    return render_to_response('supplementary/pqq.html')

def serve_raw_data(request, file_type=None, pqq_type=None, extension=None):

    # Fine the current ./data directory no matter where this call
    # is executed from (staging, production, or wherever)
    data_root = os.path.abspath("%s/data" % os.path.split(__file__)[0])

    # The URL pattern would not have matched, and thus we would not get
    # to this view if any of the keyword arguments are None
    filename = os.path.join(data_root, "%s_%s.%s" % (file_type,
                                                     pqq_type, extension))
    #from django.http import HttpResponse
    #r = HttpResponse()
    #r.write(filename)

    return serve_file(filename=filename)

    #serve_file(filename="xx")
    return r
