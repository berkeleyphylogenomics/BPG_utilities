"""Dynamic URLs for apps that can be queued"""

from django.conf.urls.defaults import *
from django.conf import settings

apps_pattern = '|'.join(settings.QUEUED_APPS)

urlpatterns = patterns('pfacts003.queued.views',
    (r'^(%s)/$' % apps_pattern, 'q_index'),
    (r'^(%s)/([\w-]+)/$' % apps_pattern, 'q_general_results'),
    (r'^(%s)/([\w-]+)/([\w-]+)/$' % apps_pattern, 'q_specific_results'),
    (r'^(%s)/([\w-]+)/([\w-]+(\.[\w-]{1,})?)$' % apps_pattern, 'q_download'),
)
