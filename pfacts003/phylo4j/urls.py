from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.phylo4j.views',
    (r'^$', 'index'),
    (r'^ppi/$', 'ppi'),
    (r'^explorer/$','explorer'),
)
