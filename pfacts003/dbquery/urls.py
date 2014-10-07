from django.conf.urls.defaults import *

urlpatterns = patterns('',
    (r'^$', 'django.views.generic.simple.direct_to_template',
          {'template': 'dbquery/index.html'}),
    (r'^help/$', 'django.views.generic.simple.direct_to_template',
          {'template': 'dbquery/help.html'}),
    (r'^about/$', 'django.views.generic.simple.direct_to_template',
          {'template': 'dbquery/about.html'}),
    (r'^tutorial/$', 'django.views.generic.simple.direct_to_template',
          {'template': 'dbquery/tutorial.html'}),
    (r'^(go|ec)/$', 'pfacts003.dbquery.views.dbquery'),
)
