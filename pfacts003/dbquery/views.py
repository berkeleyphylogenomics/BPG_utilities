from forms import GoForm, ECForm
from django.shortcuts import render_to_response
from django.core.paginator import Paginator, InvalidPage, EmptyPage

def dbquery(request, term_type):
    
    form_type = {
        'go': GoForm,
        'ec': ECForm,
    }.get(term_type)

    response_dict = {'term_type': term_type}
    if len(request.GET):
        form = form_type(request.GET)
        if form.is_valid():
            
            try:
                page = int(request.GET.get('page', '1'))
            except ValueError:
                page = 1

            ctn, pog, oog = form.results

            paginator = Paginator([
                dict(
                    canonical_tree_node=c,
                    phogs=pog[c],
                    orthologs=dict((k,list(v)) for k,v in oog[c].items())
                ) for c in ctn
            ], 20)
            response_dict['results'] = paginator.page(
                page <= 1 and 1 or \
                paginator.num_pages <= page and paginator.num_pages or \
                page
            )
                
    else:
        form = form_type()
    response_dict['form'] = form

    return render_to_response('dbquery/dbquery.html', response_dict)
