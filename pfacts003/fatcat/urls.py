from django.conf import settings
from django.conf.urls.defaults import *
from pfacts003.fatcat.views import index, results, tree, family_tree

urlpatterns = patterns('pfacts003.fatcat.views',
    (r'^family/(?P<id>.+)/tree/$', 'family_tree'),
    (r'^$', 'index'),
    (r'^(?P<id>.+)/tree/(?P<type>.+)/$', 'tree'),
    (r'^(?P<id>.+)/orthologs/$', 'orthologs_download'),
    (r'^(?P<id>.+)/alignment/$', 'alignment'),
    (r'^(?P<id>.+)/$', 'results'),
)
