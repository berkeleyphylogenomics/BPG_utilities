from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.phylobuilder.views',
    (r'^$', 'index'),
    #(r'^(\d+)/$', 'viewresults'),
)
