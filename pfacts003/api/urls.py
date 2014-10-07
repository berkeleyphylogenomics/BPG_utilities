# urls.py

from django.conf.urls.defaults import *
from piston.resource import Resource
from pfacts003.api.handlers import *
from pfacts003.protein_analysis.api import *
from pfacts003.intrepid.api import *
from pfacts003.phylo4j.api import *

sequence_classifier_handler = Resource(HMMBLASTSequenceClassifierHandler)
sequence_handler = Resource(SequenceHandler)
solr_uniprot_handler = Resource(SolrUniProtHandler)
solr_family_handler = Resource(SolrFamilyHandler)
solr_phog_handler = Resource(SolrPhogHandler)
phylofacts_search_handler = Resource(PhyloFactsSearchHandler)
fatcat_handler = Resource(FatCat)
fatcat_results_handler = Resource(FatCatResults)
fatcat_summary_handler = Resource(FatCatSummary)
fatcat_orthologs_handler = Resource(FatCatOrthologs)
fatcat2_handler = Resource(Fatcat2)
fatcat2_family_matches_handler = Resource(Fatcat2FamilyMatches)
fatcat2_functions_handler = Resource(Fatcat2Functions)
fatcat2_all_functions_handler = Resource(Fatcat2AllFunctions)
fatcat2_candidate_orthologs_handler = Resource(Fatcat2CandidateOrthologs)
fatcat2_all_candidate_orthologs_handler = Resource(Fatcat2AllCandidateOrthologs)
fatcat2_other_sequence_matches_handler = Resource(Fatcat2OtherSequenceMatches)
fatcat2_enclosing_clades_handler = Resource(Fatcat2EnclosingClades)
fatcat2_distant_clades_handler = Resource(Fatcat2DistantClades)
fatcat2_tree_handler = Resource(Fatcat2Tree)
fatcat2_go_annotation_handler = Resource(Fatcat2GOAnnotation)
fatcat2_cluster_handler = Resource(Fatcat2Cluster)
functional_site_prediction_handler = Resource(FunctionalSitePredictor)
phyloscope_ajax_request_handler = Resource(PhyloScopeAJAXRequestHandler)
phylofacts_login_handler = Resource(PhyloFactsLoginHandler)
phylofacts_logout_handler = Resource(PhyloFactsLogoutHandler)
phylofacts_user_creation_handler = Resource(PhyloFactsUserCreationHandler)
phylofacts_password_change_handler = Resource(PhyloFactsPasswordChangeHandler)
phylofacts_password_reset_handler = Resource(PhyloFactsPasswordResetHandler)

#############################################################################
# handlers for protein analysis api functions
#############################################################################
protein_analysis_job_handler = Resource(ProteinAnalysis)
protein_domain_job_handler = Resource(ProteinDomain)
# Intrepid api functions
intrepid_job_handler = Resource(Intrepid)
# for phylo4j
phylo4j_explorer_search_handler = Resource(Phylo4jExplorerSearch)
phylo4j_explorer_data_handler = Resource(Phylo4jExplorerData)

urlpatterns = patterns('',
                (r'^login/$', phylofacts_login_handler),
                (r'^logout/$', phylofacts_logout_handler),
                (r'^user_create/$', phylofacts_user_creation_handler),
                (r'^password_change/$', phylofacts_password_change_handler),
                (r'^password_reset/$', phylofacts_password_reset_handler),
                (r'^sequence_classifier_jobs/(?P<id>(\w+))/$', sequence_classifier_handler),
                (r'^sequence_classifier_jobs/$', sequence_classifier_handler),
                (r'^sequence/$', sequence_handler),
                (r'^solr_uniprot/$', solr_uniprot_handler),
                (r'^solr_family/$', solr_family_handler),
                (r'^solr_phog/$', solr_phog_handler),
                (r'^phylofacts/search/$', phylofacts_search_handler),
                # For Phyloscope ajax request
                (r'^phyloscope/$', phyloscope_ajax_request_handler),
                # Functional Site Prediction Jobs
                (r'^functional_site_prediction_jobs/(?P<jobID>(\d+))/$', functional_site_prediction_handler),
                (r'^functional_site_prediction_jobs/$', functional_site_prediction_handler),
                #Obese Kitten
                (r'^fatcat/$', fatcat_handler),
                (r'^fatcat/(?P<fatcat_id>(\d+))/$', fatcat_handler),
                (r'^fatcat/(?P<fatcat_id>(\d+))/summary/$', fatcat_summary_handler),
                (r'^fatcat/(?P<fatcat_id>(\d+))/results/$', fatcat_results_handler),
                (r'^fatcat/(?P<fatcat_id>(\d+))/orthologs/$', fatcat_orthologs_handler),
                # new fatcat
                (r'^fatcat2/$', fatcat2_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/$', fatcat2_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/all_candidate_orthologs/$', fatcat2_all_candidate_orthologs_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/all_functions/$', fatcat2_all_functions_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/functions/$', fatcat2_functions_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/candidate_orthologs/$', fatcat2_candidate_orthologs_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/other_sequence_matches/$', fatcat2_other_sequence_matches_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/family_matches/$', fatcat2_family_matches_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/enclosing_clades/$', fatcat2_enclosing_clades_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/distant_clades/$', fatcat2_distant_clades_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/tree/$', fatcat2_tree_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/annotations/go/(?P<annotation_id>(.+))/', fatcat2_go_annotation_handler),
                (r'^fatcat2/(?P<fatcat_id>(\d+))/cluster/(?P<cluster_id>(\d+))/', fatcat2_cluster_handler),
                
                (r'^SAS/$', protein_analysis_job_handler),
                (r'^SAS/(?P<id>(\d+))/$', protein_analysis_job_handler),
                (r'^SAS/domain/(?P<id>(\d+))/$', protein_domain_job_handler),
                (r'^intrepid/$', intrepid_job_handler),
                (r'^intrepid/(?P<id>(\d+))/$', intrepid_job_handler),
                # for phylo4j api urls
                # phylo4j explorer
                (r'phylo4j/explorer/search/$', phylo4j_explorer_search_handler),
                (r'phylo4j/explorer/(?P<id>(\d+))/(?P<t>(.+))/$', phylo4j_explorer_data_handler),
) 
