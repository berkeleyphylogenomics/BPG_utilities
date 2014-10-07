from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.common.views',
    (r'^clean_webserver_temp$', 'clean_webserver_temp'),
)
