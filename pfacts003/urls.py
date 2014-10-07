from django.conf import settings
from django.conf.urls.defaults import *
from pfacts003.views import member, members, kimmen_personal, temp_out_of_service

urlpatterns = patterns('',
    (r'^$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/index_new.html'}),
    (r'^publications/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/publications.html'}),
    (r'^software/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/software.html'}),
    (r'^positions/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/positions.html'}),                       
    (r'^research/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/research.html'}),
    (r'^affiliations/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/affiliations.html'}),
    (r'^funding/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/funding.html'}),
    (r'^homology_model/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/homology_model.html'}),

    (r'^contact_us/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'common/contact_us.html'}),
    (r'^courses/$', 'django.views.generic.simple.direct_to_template',
        {'template': 'home/courses.html'}),
    (r'^members/$', members),
    (r'^members/(?P<member>\w+)/$', member),
    (r'^INTREPID/', include('pfacts003.intrepid.urls')),
    (r'^intrepid/', include('pfacts003.intrepid.urls')),                       
    #(r'^FSP/', include('pfacts003.FSP.urls')),
    (r'^SCI-PHY/', temp_out_of_service),
    (r'^flowerpower/', temp_out_of_service),
    (r'^phylobuilder/', temp_out_of_service),
    (r'^(?i)phylobuild/', temp_out_of_service),                       
    (r'^members/kimmen/(?P<page>\w+\.(html|pdf))/$', kimmen_personal),
                       
    # Legacy support redirects
    # These should never be deprecated as they are urls that
    # have been published in scientific journals
    # NAR 2010 1-6 doi:10.1093/nar/gkq298
    (r'^satchmo(?:-js|2)/$', 'django.views.generic.simple.redirect_to', 
                             {'url': '/q/satchmo/'}),
    (r'^SATCHMO(?:-JS|2)/$', 'django.views.generic.simple.redirect_to', 
                             {'url': '/q/satchmo/'}),                       
    # NAR 2010 1-6 doi:10.1093/nar/gkq298
    (r'^satchmo-js/supplementary/',
        'django.views.generic.simple.redirect_to', 
             {'url': '/supplementary/satchmo/'}),

    # Supportive and shared apps
    (r'^q/', include('pfacts003.queued.urls')),
    (r'^supplementary/', include('pfacts003.supplementary.urls')),

    # Staff only applications
    (r'^staff/', include('pfacts003.common.staff_urls')),

    # Each app should have it's own urls so things are more modular and 
    # configurable. These in production
    (r'^home/', include('pfacts003.home.urls')),
    (r'^(?i)satchmo/', include('pfacts003.satchmo.urls')),
    #(r'^phog/', include('pfacts003.phog.urls')),
    (r'^phog/', 'django.views.generic.simple.redirect_to', {'url': '/fatcat/'}),
    (r'^sciphy/', include('pfacts003.sciphy.urls')),
    (r'^phylofacts/', include('pfacts003.phylofacts.urls')),
    (r'^phyloscope/', include('pfacts003.phyloscope.urls')),
    (r'^api/', include('pfacts003.api.urls')),
    (r'^HYPNO/', include('pfacts003.hypno.urls')),
    #(r'^SAS/', include('pfacts003.protein_analysis.urls')), 
    (r'^fatcat/', include('pfacts003.fatcat.urls')),
    (r'^phylo4j/', include('pfacts003.phylo4j.urls')),
)


# The following is for development use only. Please always ensure
# DEBUG is  False in production (as it should be by default).

if settings.DEBUG is True:
    urlpatterns += patterns('',
        # Each app should have it's own urls so things are more modular and 
        # configurable. These in test
        (r'^dbquery/', include('pfacts003.dbquery.urls')),
        (r'^phylobuilder/', include('pfacts003.phylobuilder.urls')),
        (r'^kegg/', include('pfacts003.kegg.urls')),

    )
    #urlspatterns += staticfiles_urlpatterns()

if settings.DEBUG is True:
    urlpatterns += patterns('',
        # For static serving on local debug server
        url(r'^static/(?P<path>.*)$', 'django.views.static.serve',
            {'document_root': settings.MEDIA_ROOT,
             'show_indexes': True}
        ),
    )
       
