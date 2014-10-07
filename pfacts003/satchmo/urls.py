from django.conf.urls.defaults import *

urlpatterns = patterns('',
    (r'^$', 'pfacts003.satchmo.views.index'),
    (r'^about/$', 'pfacts003.satchmo.views.about'),
    (r'^help/$', 'pfacts003.satchmo.views.help'),
)
