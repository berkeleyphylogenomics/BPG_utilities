from django.conf.urls.defaults import *

_file_types = "homologs|models|MSA"
_file_extensions = "txt|pdb|afa"

urlpatterns = patterns('pfacts003.supplementary.views',
    (r'^pqq/$', 'index_pqq'),

    # Accept patterns like:
    # .....pqq/homologs_pqqA.txt
    # .....pqq/models_pqqB.pdb
    (r'^pqq/(?P<file_type>%s)_(?P<pqq_type>pqq[A-F])\.(?P<extension>%s)$' % (
        _file_types, _file_extensions), 'serve_raw_data'),
)


urlpatterns += patterns('',
    (r'satchmo/$',
          'django.views.generic.simple.direct_to_template',
              {'template': 'supplementary/satchmo/supplementary_material.html'}),
    (r'satchmo/webserver/$',
          'django.views.generic.simple.direct_to_template',
            {'template': 'supplementary/satchmo/supplementary_material_webserver.html'}),
    (r'satchmo/algorithm/$',
          'django.views.generic.simple.direct_to_template',
            {'template': 'supplementary/satchmo/supplementary_material_algorithm.html'}),
)
