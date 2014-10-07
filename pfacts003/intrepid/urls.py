from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.intrepid.views',
    (r'^$', 'index'),
    (r'^pdb/(?P<pid>.+)/$', 'serve_pdb'),
    (r'^(?P<id>.+)/download/(?P<type>.+)/$', 'download'),
    (r'^(?P<id>.+)/pdb/$', 'serve_this_job_pdb'),
    (r'^(?P<id>.+)/$', 'results'),
)
