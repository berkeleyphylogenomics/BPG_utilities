from django.conf import settings
from django.conf.urls.defaults import *
from pfacts003.phog.views import alignment_blast_pairwise, \
                              approximate_matches, \
                              orthoscope, orthoscope_get, \
                              orthologs, orthologs_as_csv, \
                              orthologs_quickstart, orthologs_faq, \
                              orthologs_about, orthologs_assessment, \
                              orthologs_supplement, orthologs_tutorial, \
                              phog, phog_as_csv, orthologs_downloads, \
                              phog_tree, phog2phog, hyperscope, netscope, \
                              interactors, interactors_as_csv, phog_tree_json
from pfacts003.utils.id_patterns import uniprot_accession_pat1, uniprot_accession_pat2, \
  uniprot_taxon_pat, swissprot_desc_pat, scopid_pat, bpgid_pat, phog_pat, \
  gi_pat, leftid_pat, phog_tight_pat, phog_medium_pat, phog_loose_pat, \
  phog_custom_pat1, phog_custom_pat2, phog_custom_pat3, phog_custom_pat4, \
  taxon_id_pat

urlpatterns = patterns('',
    # These patterns are already being referenced by in publications
    # and other linked pages. They need to stay in production (but can
    # be refactored).
    (r'^orthologs/supplement/(v1)/$', orthologs_supplement),
    (r'^orthologs/downloads/', orthologs_downloads),
    (r'^orthologs/quickstart/$', orthologs_quickstart),
    (r'^orthologs/faq/$', orthologs_faq),
    (r'^orthologs/FAQ/$', orthologs_faq),
    (r'^orthologs/tutorial/$', orthologs_tutorial),
    (r'^orthologs/assessment/$', orthologs_assessment),
    (r'^alignment/$', alignment_blast_pairwise),
    (r'^orthoscope/(%sd?)/$' % bpgid_pat, orthoscope),
    (r'^orthoscope/?$', orthoscope_get),
    (r'^orthologs/$', orthologs),
    (r'^orthologs/about/$', orthologs_about),
    (r'^orthologs/(%s)/$' % uniprot_accession_pat1, orthologs),
    (r'^orthologs/(%s_%s)/$' % (uniprot_accession_pat1, uniprot_taxon_pat),
      orthologs),
    (r'^orthologs/(%s)/$' % uniprot_accession_pat2, orthologs),
    (r'^orthologs/(%s_%s)/$' % (uniprot_accession_pat2, uniprot_taxon_pat),
      orthologs),
    (r'^orthologs/(%s_%s)/$' % (swissprot_desc_pat, uniprot_taxon_pat),
      orthologs),
    (r'^orthologs/(%s)/$' % gi_pat, orthologs),
    (r'^orthologs/csv/(%s)/$' % uniprot_accession_pat1, orthologs_as_csv),
    (r'^orthologs/csv/(%s_%s)/$' % (uniprot_accession_pat1, 
                                    uniprot_taxon_pat), orthologs_as_csv),
    (r'^orthologs/csv/(%s)/$' % uniprot_accession_pat2, orthologs_as_csv),
    (r'^orthologs/csv/(%s_%s)/$' % (uniprot_accession_pat2, 
                                    uniprot_taxon_pat), orthologs_as_csv),
    (r'^orthologs/csv/(%s_%s)/$' % (swissprot_desc_pat, uniprot_taxon_pat),
      orthologs_as_csv),
    (r'^orthologs/csv/(%s)/$' % gi_pat, orthologs_as_csv),

    (r'^$', orthologs),
    (r'^(%s)/$' % phog_pat, phog),
    (r'^csv/(%s)/$' % phog_pat, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_pat, phog_pat), phog2phog),
    (r'^hyper/(%s)/$' % phog_pat, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_pat, taxon_id_pat), hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_pat, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_pat, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_pat, phog_tree),
    (r'^tree/(%s).json/$' % phog_pat, phog_tree_json),
    (r'^(%s)/$' % phog_tight_pat, phog),
    (r'^csv/(%s)/$' % phog_tight_pat, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_tight_pat, phog_tight_pat), phog2phog),
    (r'^hyper/(%s)/$' % phog_tight_pat, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_tight_pat, taxon_id_pat), hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_tight_pat, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_tight_pat, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_tight_pat, phog_tree),
    (r'^(%s)/$' % phog_medium_pat, phog),
    (r'^csv/(%s)/$' % phog_medium_pat, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_medium_pat, phog_medium_pat), 
                              phog2phog),
    (r'^hyper/(%s)/$' % phog_medium_pat, hyperscope),
    (r'^hyper/(%s)_(%s)/$'% (phog_medium_pat, taxon_id_pat), hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_medium_pat, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_medium_pat, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_medium_pat, phog_tree),
    (r'^(%s)/$' % phog_loose_pat, phog),
    (r'^csv/(%s)/$' % phog_loose_pat, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_loose_pat, phog_loose_pat), phog2phog),
    (r'^hyper/(%s)/$' % phog_loose_pat, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_loose_pat, taxon_id_pat), hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_loose_pat, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_loose_pat, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_loose_pat, phog_tree),
    (r'^(%s)/$' % phog_custom_pat1, phog),
    (r'^csv/(%s)/$' % phog_custom_pat1, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_custom_pat1, phog_custom_pat1), 
                                  phog2phog),
    (r'^hyper/(%s)/$' % phog_custom_pat1, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_custom_pat1, taxon_id_pat), 
                                hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_custom_pat1, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_custom_pat1, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_custom_pat1, phog_tree),
    (r'^(%s)/$' % phog_custom_pat2, phog),
    (r'^csv/(%s)/$' % phog_custom_pat2, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_custom_pat2, phog_custom_pat2), 
                                  phog2phog),
    (r'^hyper/(%s)/$' % phog_custom_pat2, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_custom_pat2, taxon_id_pat), 
                                hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_custom_pat2, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_custom_pat2, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_custom_pat2, phog_tree),
    (r'^(%s)/$' % phog_custom_pat3, phog),
    (r'^csv/(%s)/$' % phog_custom_pat3, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_custom_pat3, phog_custom_pat3), 
                                  phog2phog),
    (r'^hyper/(%s)/$' % phog_custom_pat3, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_custom_pat3, taxon_id_pat), 
                                hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_custom_pat3, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_custom_pat3, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_custom_pat3, phog_tree),
    (r'^(%s)/$' % phog_custom_pat4, phog),
    (r'^csv/(%s)/$' % phog_custom_pat4, phog_as_csv),
    (r'^phog/(%s)_(%s)/$' % (phog_custom_pat4, phog_custom_pat4), 
                                  phog2phog),
    (r'^hyper/(%s)/$' % phog_custom_pat4, hyperscope),
    (r'^hyper/(%s)_(%s)/$' % (phog_custom_pat4, taxon_id_pat), 
                                hyperscope),
    (r'^tree/(%s)_(%s)_(%s)/$' % (phog_custom_pat4, leftid_pat, leftid_pat), 
      phog_tree),
    (r'^tree/(%s)_(%s)/$' % (phog_custom_pat4, leftid_pat), phog_tree),
    (r'^tree/(%s)/$' % phog_custom_pat4, phog_tree),
    
    (r'^ppi/$', interactors),
    (r'^ppi/(%s)/$' % uniprot_accession_pat1, interactors),
    (r'^ppi/(%s_%s)/$' % (uniprot_accession_pat1, uniprot_taxon_pat),
      interactors),
    (r'^ppi/(%s)/$' % uniprot_accession_pat2, interactors),
    (r'^ppi/(%s_%s)/$' % (uniprot_accession_pat2, uniprot_taxon_pat),
      interactors),
    (r'^ppi/(%s_%s)/$' % (swissprot_desc_pat, uniprot_taxon_pat),
      interactors),
    (r'^ppi/(%s)/$' % gi_pat, interactors),
    (r'^ppi/csv/(%s)/$' % uniprot_accession_pat1, interactors_as_csv),
    (r'^ppi/csv/(%s_%s)/$' % (uniprot_accession_pat1, 
                                    uniprot_taxon_pat), interactors_as_csv),
    (r'^ppi/csv/(%s)/$' % uniprot_accession_pat2, interactors_as_csv),
    (r'^ppi/csv/(%s_%s)/$' % (uniprot_accession_pat2, 
                                    uniprot_taxon_pat), interactors_as_csv),
    (r'^ppi/csv/(%s_%s)/$' % (swissprot_desc_pat, uniprot_taxon_pat),
      interactors_as_csv),
    (r'^ppi/csv/(%s)/$' % gi_pat, interactors_as_csv),

    (r'^net/$', netscope),

    )

