from django.shortcuts import render_to_response

def research(request):
    return render_to_response("home/research.html",
            {'main_viewing': 'home', 'sub_viewing': 'research'})

def publications(request):
    return render_to_response("home/publications.html",
            {'main_viewing': 'home', 'sub_viewing': 'publications'})

def affiliations(request):
    return render_to_response("home/affiliations.html",
            {'main_viewing': 'home', 'sub_viewing': 'affiliations'})

def members(request):
    return render_to_response("home/members.html",
            {'main_viewing': 'home', 'sub_viewing': 'members'})

def seminars(request):
    return render_to_response("home/seminars.html",
            {'main_viewing': 'home', 'sub_viewing': 'seminars'})

def maintenance(request):
    return render_to_response("maintenance.html")
