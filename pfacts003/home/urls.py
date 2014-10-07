from django.conf.urls.defaults import *

urlpatterns = patterns('pfacts003.home.views',
    (r'^research', 'research'),
    (r'^publications', 'publications'),
    (r'^affiliations', 'affiliations'),
    (r'^members', 'members'),
    (r'^seminars', 'seminars'),
    (r'', 'research'),
    (r'^$', 'research'),
)
