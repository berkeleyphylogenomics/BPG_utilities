from django.shortcuts import render_to_response
from django.template import RequestContext

def members(request):
    return render_to_response('home/members.html', {}, context_instance=RequestContext(request))

def member(request, member):
    template = "home/members/%s.html" % member
    return render_to_response(template)

def kimmen_personal(request, page):
    template = "home/members/kimmen_personal/%s" % page
    return render_to_response(template)

def temp_out_of_service(request):
    template = "home/temp_out_of_service.html"
    return render_to_response(template)
