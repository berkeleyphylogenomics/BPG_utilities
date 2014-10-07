from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.protein_analysis.views',
    (r'^$', 'index'),
    (r'^(?P<id>.+)/$', 'results'),
)
