from django.conf import settings
from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.hypno.views',
    (r'^$', 'home'),
    (r'^data/source/$', 'serve_hypno_source_code'),
    (r'^data/supplementary/$', 'serve_hypno_supplementary_data'),
    (r'^help/$', 'help'),
)
