from django.conf.urls.defaults import *
from pfacts003.phyloscope.views import phyloscope, phyloscope_get, \
                                        phyloscope_quick_start
from pfacts003.utils.id_patterns import bpgid_pat

urlpatterns = patterns('',
    (r'^\?$', phyloscope_get),
    (r'^(%s)/$' % bpgid_pat, phyloscope),
    (r'^help/$', phyloscope_quick_start),
    (r'^quickstart/$', phyloscope_quick_start),
)
