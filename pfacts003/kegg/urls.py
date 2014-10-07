from django.conf.urls.defaults import *
from pfacts003.kegg.views import kegg_map
from pfacts003.utils.id_patterns import kegg_map_pat, taxon_id_pat

urlpatterns = patterns('',
    (r'^(%s)/$' % kegg_map_pat, kegg_map),
    (r'^(%s)/taxon(%s)/$' % (kegg_map_pat, taxon_id_pat), kegg_map),

    # Uncomment this for admin:
#     (r'^admin/', include('django.contrib.admin.urls')),
)
