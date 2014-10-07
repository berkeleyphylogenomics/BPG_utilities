import os, string
from django.utils.safestring import mark_safe
from id_patterns import bpgid_pat, scopid_pat
from utils import remove_duplicates  
from urllib import quote_plus

"""
  Tooltip functions and tags ("tooltip_val") are defined in a JavaScript file elsewhere, 
    for instance, /media/js/phog_tooltip.js
  Icons are image files in /media/icons/.  The class "icon" removes the border and centers 
    the icon image.
"""

# Generate the HTML JavaScript calls to show and hide the tooltip window.
# The functions are in phog_tooltip.js and bpg_tooltip.js.
# Tooltip display function is showToolTip by default; showToolTip requires a tag that
#   specifies which tooltip text to display.
def make_tooltip(tooltip_val, tooltip_func):
  if tooltip_val:
    if not tooltip_func:
      tooltip_func = "showToolTip";
    return """onmouseover="%s('%s')" onmouseout="UnTip()" """ % (tooltip_func, str(tooltip_val))
  else:
    return ''

# Generate the HTML for a text link to a URL with an optional tooltip tag.
def make_link(text, url, tooltip_val = "", tooltip_func = "", class_str = ""):
  if text and url:
    if class_str:
      class_str = 'class="%s" ' % class_str
    return mark_safe('<a %s%shref="%s">%s</a>' % (class_str, make_tooltip(tooltip_val, tooltip_func), url, text))
  else:
    return ''

# Generate tooltipped text without a link to a URL as used for column headings, etc.  
def make_tooltipped_text(text, tooltip_val, tooltip_func = "", class_str = ""):
  if class_str:
    class_str = "class='%s' " % class_str
  if text:
    return mark_safe('<span %s%s>%s</span>' % (class_str, make_tooltip(tooltip_val, tooltip_func), text))
  else:
    return ''

# Generates the HTML to display an icon from the Django media files, /media/icons/.  
# Class "icon" suppresses border and centers icon.
def make_icon(icon):
  return mark_safe("""<img class="icon" src="/static/img/icons/%s" />""" \
          % icon)

# Generate the HTML for an icon linked to a URL with optional tooltip tag.
def make_icon_link(icon, url, tooltip_val = "", tooltip_func = ""):
  if url: 
    return make_link(make_icon(icon), url, tooltip_val, tooltip_func)
  else:
    return ''
  
"""
  Definitions of URLs and links follow in alphabetical order...
"""  

def alignment_blast_pairwise_url(identifier, file):
  return ('/phog/alignment/?identifier=%s&file=%s&location=phogblast' % (identifier, file))

def alignment_blast_pairwise_link(text, identifier, file):
  return make_link(text, alignment_blast_pairwise_url(identifier, file), 'blast_alignment')

def alignment_hmm_pairwise_url(identifier, basename, query_title, query_length):
  return('/phylofacts/alignment/?identifier=%s&basename=%s&qdesc=%s&qlen=%d') \
      % (identifier, basename, quote_plus(query_title), query_length)

def alignment_hmm_pairwise_link(text, identifier, query_title, query_length,
                                  basename):
  return make_link(text, alignment_hmm_pairwise_url(identifier, basename,
                                          query_title, query_length),
                    'blast_alignment')

def biocyc_url(id):
  if id:
    return 'http://biocyc.org/getid?id=%s' % id
  else:
    return ''

def biocyc_icon(id):
  return make_icon_link('icon_biocyc_12.png', biocyc_url(id), 'biocyc')

def family_url(family):  
  if isinstance(family, basestring):
      id_string = family
  elif isinstance(family, int):
      id_string = "bpg%07d" % family
  else:
      id_string = "bpg%07d" % family.id
      
  return ('/phylofacts/family/bpg%s/' % id_string)

def family_link(family):
  family = str(family)
  return make_link(family, family_url(family), 'family')

def family_orthoscope_url(family):  
  family = str(family)
  if family:
    return ('/phog/orthoscope/?family=%s' % family)
  else:
    return ''

def family_orthoscope_link(family):
  family = str(family)
  return make_link(family, family_orthoscope_url(family), 'phyloscope') 

def family_orthoscope_link_in_phog_table(family, n_seqs, text = 'Full'):
  family = str(family)
  if n_seqs < 300:
    return make_link(text, family_orthoscope_url(family), n_seqs, 'showNSeqsTip') 
  else:
    return ''

def family_phyloscope_url(family):  
  family = str(family)
  if family:
    return ('/phyloscope/%s/' % family)
  else:
    return ''

def family_phyloscope_link(family):
  family = str(family)
  return make_link(family, family_phyloscope_url(family), 'phyloscope') 

def dip_url(dip_id):
  if dip_id:
    return 'http://dip.doe-mbi.ucla.edu/dip/DIPview.cgi?ID=%s' % dip_id
  else:
    return ''

def dip_link(dip_id):
  return make_link(dip_id, dip_url(dip_id))

def kegg_url(id):
  if id:
    return 'http://www.genome.jp/dbget-bin/www_bget?%s' % id
  else:
    return ''

def kegg_icon(id):
  return make_icon_link('icon_japan_flag_18.png', kegg_url(id), 'kegg')

def ncbi_entrez_url(identifier):
  return 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=%s' % identifier

def ncbi_entrez_link(identifier):
  return make_link(identifier, ncbi_entrez_url(identifier), 'ncbi') 

def ncbi_taxonomy_url(ncbi_tax_id):
  return 'http://www.uniprot.org/taxonomy/%s' % ncbi_tax_id

def ncbi_taxonomy_link(ncbi_tax, ncbi_tax_id):
  species = str(ncbi_tax).replace(" (", "<br />(")
  return make_link(species, ncbi_taxonomy_url(ncbi_tax_id), 'species')

def netscope_link(uniprot_id, ortholog_type = 0, threshold = 0.0, text = 'PPI'):
  return make_link(text, netscope_url(uniprot_id, ortholog_type, threshold), 'ppi')

def netscope_url(id, ortholog_type, threshold):
  if id:
    return '/phog/net/?sequence_id=%s&ortholog_type=%s&threshold=%s' % (id, ortholog_type, threshold)
  else:
    return ''

def phylofacts_pfam_url(pfam_acc):
    pfam_acc = str(pfam_acc)
    return '/phylofacts/pfam/%s' % pfam_acc

def pfam_url(pfam):
  pfam = str(pfam)
  return 'http://pfam.sanger.ac.uk/family?id=%s' % pfam

def pfam_link(pfam_name, pfam_description, range=None):
#  if range:
#    start, end = range
#    anchor = '%s (%d, %d)' % (pfam_description, start, end)
#  else:
#    anchor=pfam_description
  anchor = pfam_name
  return make_link(anchor, pfam_url(pfam_name), 'pfam')

def pf_pfam_link(pfam_name, pfam_acc):
    anchor = pfam_name
    return make_link(anchor, phylofacts_pfam_url(pfam_acc), 'pfam')

def pmid_url(pmid):
  return "http://ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=search&term=%d" % pmid

def pmid_link(pmid):
  return make_link('PMID:%d' % pmid, pmid_url(pmid), 'PubMed')

def kegg_map_url(kegg_map):
  id_str = "map%05d" % kegg_map.id
  return 'http://www.genome.jp/kegg/pathway/map/%s.html' % id_str

def kegg_map_link(kegg_map):
  return make_link(kegg_map.title, kegg_map_url(kegg_map), 'kegg')

def phog_tree_url(phog_accession, ortholog_leftid = "", query_leftid = ""):
  ortholog_str = ""
  query_str = ""
  if ortholog_leftid:
    ortholog_str = "_%s" % str(ortholog_leftid)
  if query_leftid:
    query_str = "_%s" % str(query_leftid)
  return '/phog/tree/%s%s%s/' % (phog_accession, ortholog_str, query_str)

def phog_tree_link(phog_accession, ortholog_leftid = "", query_leftid = "", text = 'PHOG'):
  return make_link(text, phog_tree_url(phog_accession, ortholog_leftid, query_leftid), 'phog_tree')  

def sequence_orthologs_url(seq):
  seq = str(seq)
  if seq:
    return '/phog/orthologs/%s/' % seq
  else:
    return ''

def sequence_orthologs_link(seq, text = 'Orthologs'):
  return make_link(text, sequence_orthologs_url(seq), 'orthologs') 

def swissprot_icon(id):
  return make_icon_link('icon_swiss_flag_12.png', uniprot_url(id), 'in_swissprot')

def uniprot_url(uniprot):
  uniprot = str(uniprot)
  return 'http://www.uniprot.org/uniprot/%s' % uniprot

def uniprot_link(uniprot):
  uniprot = str(uniprot)
  return make_link(uniprot, uniprot_url(uniprot), 'gene')

def uniprot_ec_url(ec):
  if ec:
    e_class = string.replace(str(ec),".-", "")
    return 'http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=%s' % e_class
  else:
    return ''

def uniprot_ec_link(ec, experimental=False):
  if experimental:
    return make_link(ec, uniprot_ec_url(ec), 'ec') + mark_safe(' * ')
  else:
    return make_link(ec, uniprot_ec_url(ec), 'ec')

def uniprot_evidence_url(uniprot):
  uniprot = str(uniprot)
  return 'http://www.uniprot.org/uniprot/%s#section_terms' % uniprot

def uniprot_evidence_icon(uniprot):
  uniprot = str(uniprot)
  return make_icon_link('icon_flask_12.png', uniprot_evidence_url(uniprot), 'go_evidence')

def uniprot_lit_url(uniprot):
  uniprot = str(uniprot)
  return 'http://www.uniprot.org/uniprot/%s#section_ref' % uniprot

def uniprot_lit_icon(uniprot):
  uniprot = str(uniprot)
  return make_icon_link('icon_literature_12.png', uniprot_lit_url(uniprot), 'literature')

def uniprot_orthologs_url(uniprot):
  uniprot = str(uniprot)
  return '/phog/orthologs/%s/' % uniprot

def uniprot_orthologs_link(uniprot, text = 'Orthologs'):
  uniprot = str(uniprot)
  return make_link(text, uniprot_orthologs_url(uniprot), 'orthologs') 

def uniprot_ppi_link(uniprot_id, text = 'PPI'):
  return make_link(text, uniprot_ppi_url(uniprot_id), 'ppi')

def uniprot_ppi_url(id):
  if id:
    return '/ppi/%s/' % id
  else:
    return ''

"""
  Functions to parse Pfam designations from pfam_color_file and create links
"""

def make_kegg_map_links(kegg_maps):
  kegg_map_links = [kegg_map_link(kegg_map) for kegg_map in kegg_maps]
  return mark_safe(', '.join(kegg_map_links))

def make_uniprot_pfam_links(uniprot):
  pfams = set([pfama_pfamseq.pfamA for pfama_pfamseq 
                in uniprot.get_pfam_hits()])
  pfam_links = [pf_pfam_link(pfam.pfama_id, pfam.pfama_acc)
                          for pfam in pfams]
  return mark_safe(', '.join(pfam_links))

def make_pfam_links(family, include_nonlinks = False):
  pfams = [(start, end, pfam) for (pfam, start, end) in family.get_pfams()]
  pfams.sort()
  included_pfams = set()
  nonredundant_pfams = []
  for start, end, pfam in pfams:
    if pfam not in included_pfams:
      nonredundant_pfams = nonredundant_pfams + [pfam]
      included_pfams.add(pfam)
  pfam_links = [pf_pfam_link(pfam.name, pfam.accession) 
                  for pfam in nonredundant_pfams]
  # TODO: 2010/06/14 If include_nonlinks, include transmembrane and signal
  # peptide predictions
  return mark_safe(' | '.join(pfam_links))

def make_ec_links(family):
    experimental_ecs = family.canonical_root_node().get_ecs(experimental=True, 
                          ortholog_type = OrthologTypes.PHOG_T_Custom,
                          threshold=10000.0)
    nonexperimental_ecs = family.canonical_root_node().get_ecs(
              experimental=False, ortholog_type = OrthologTYpes.PHOG_T_Custom,
              threshold=10000.0)
    ec_links = mark_safe(','.join([uniprot_ec_link(ec, experimental=True) 
                                      for ec in experimental_ecs] + 
                                      [uniprot_ec_link(ec, experimental=False)
                                      for ec in nonexperimental_ecs]))
    return ec_links

def make_go_summary_line(go_summary):
  return ', '.join(['<a href="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s">%s</a> (<a href="http://www.geneontology.org/GO.evidence.shtml#%s" title="%s">%s</a>)' \
    % (acc, name, evidence, description, evidence) 
    for (acc, name, evidence, description) in go_summary])

def make_go_summary_divs(go_component, long_summary, short_summary):
  result = ""
  if len(long_summary) > 0:
    if len(long_summary) != len(short_summary):
      result =  '<form name="long_%s_form">' % go_component
      result += '<div id="long_%s" class="go_summary-visible">' % go_component
      result += mark_safe(make_go_summary_line(long_summary))
      result += '<input type="button" value="(hide)" class="fake_btn"'
      result += """onclick="switchVisible('long_%s', 'short_%s')">""" \
                % (go_component, go_component)
      result += "</div>"
      result += '<div id="short_%s" class="go_summary-invisible">' \
                % go_component
      result += mark_safe(make_go_summary_line(short_summary))
      result += '<input type="button" value="(more...)" class="fake_btn"'
      result += """onclick="switchVisible('short_%s', 'long_%s')">""" \
                % (go_component, go_component)
      result += '</div>'
      result += '</form>'
      result += '<script language="JavaScript">'
      result += 'switchVisible("long_%s", "short_%s")' % (go_component, 
                                                        go_component)
      result += '</script>'
    else:
      result = mark_safe(make_go_summary_line(long_summary))
  return mark_safe(result)
