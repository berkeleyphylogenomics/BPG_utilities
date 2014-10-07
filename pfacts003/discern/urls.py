from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.discern.views',
    (r'(\d+)/$', 'viewresults'),
    (r'(\d+)/(.{5})/$', 'viewresults'),
    (r'scores/(\d+)/$', 'viewscores'),
    (r'^$', 'index'),
)
