from django.conf.urls.defaults import *

# The SCI-PHY alogorithm webserver has not yet been ported from Limahuli. When
# it is ported, this is the Django app that it will occupy. The app is
# currently here so that we can redirect the supplementary materials to
# Limahuli per the link in the /publications/ document.

urlpatterns = patterns('',
    (r'^$', 'django.views.generic.simple.redirect_to',
                    {'url': 'http://phylogenomics.berkeley.edu/SCI-PHY/supplemental/index.html'}),
    (r'^supplemental/$', 'django.views.generic.simple.redirect_to',
                    {'url': 'http://phylogenomics.berkeley.edu/SCI-PHY/supplemental/index.html'}),
)
