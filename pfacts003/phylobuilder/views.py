#from pfacts003.shared_views import doIndex, doViewresults
#from pfacts003.shared_views import doIndex
from django.shortcuts import render_to_response

def index(request):
    """handles displaying the initial phylobuilder form,
 launching new analysis jobs, and displaying
 the progress of running jobs
"""

    #return doIndex(request, 'phylobuilder', 'phylobuilder.py -i')

    return render_to_response('phylobuilder/index.html')


# Currently not importing doViewresults == removed until this can be fixed
#def viewresults(request, familyid):
#    """pulls up a display of existing results"""
#
#    return doViewresults(request, 'phylobuilder', familyid)
