from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.hpy.views',
    (r'^$', 'index'),
    (r'^about/$', 'about'),
    (r'^sequence/$', 'sequence'),
)
