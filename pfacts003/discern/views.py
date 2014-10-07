from pfacts003.phylofacts.models import *
from pfacts003.common.shared_views import  doIndex, \
                                     viewStructResults, \
                                     viewSequenceScores

from django.shortcuts import render_to_response

#
# handles displaying the initial intrepid form,
# launching new analysis jobs, and displaying
# the progress of running jobs
#

def index(request):
    return doIndex(request, 'discern', 'discern.py')


#
# pulls up a display of existing results
#

def viewresults(request, familyid, pdbid=''):
    return viewStructResults(request, familyid, pdbid, 'discern')

#
# views discern scores of existing results
#

def viewscores(request, familyid):
    return viewSequenceScores(request, familyid, 'discern')

