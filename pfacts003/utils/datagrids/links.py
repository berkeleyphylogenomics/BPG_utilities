import os
import string

from django.utils.safestring import mark_safe

from pfacts003.utils.files import open_pfam_color_file
from pfacts003.utils.id_patterns import bpgid_pat, scopid_pat
from pfacts003.utils.misc import remove_duplicates  

"""
Tooltip functions and tags ("tooltip_val") are defined in a JavaScript file
elsewhere, for instance, /media/js/phog_tooltip.js Icons are image files in
/media/icons/.  The class "icon" removes the border and centers the icon image.
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
  return ('/alignment/?identifier=%s&file=%s&location=phogblast' % (identifier, file))

def alignment_blast_pairwise_link(text, identifier, file):
  return make_link(text, alignment_blast_pairwise_url(identifier, file), 'blast_alignment')

def biocyc_url(id):
  if id:
    return 'http://biocyc.org/getid?id=%s' % id
  else:
    return ''

def biocyc_icon(id):
  return make_icon_link('icon_biocyc_12.png', biocyc_url(id), 'biocyc')

def family_url(family):  
  family = str(family)
  return ('/family/%s/' % family)

def family_link(family):
  family = str(family)
  return make_link(family, family_url(family), 'family')

def family_orthoscope_url(family):  
  family = str(family)
  if family:
    return ('/orthoscope/?family=%s' % family)
  else:
    return ''

def family_orthoscope_link(family):
  family = str(family)
  return make_link(family, family_orthoscope_url(family), 'phyloscope') 

def family_orthoscope_link_in_phog_table(family, n_seqs, text = 'Full'):
  family = str(family)
  return make_link(text, family_orthoscope_url(family), n_seqs, 'showNSeqsTip') 

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

def pfam_url(pfam):
  pfam = str(pfam)
  return 'http://pfam.sanger.ac.uk/family?id=%s' % pfam

def pfam_link(pfam):
  pfam = str(pfam)
  return make_link(pfam, pfam_url(pfam), 'pfam')

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
    return '/orthologs/%s/' % seq
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
    e_class = string.replace(ec,".-", "")
    return 'http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=%s' % e_class
  else:
    return ''

def uniprot_ec_link(ec):
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
  return '/orthologs/%s/' % uniprot

def uniprot_orthologs_link(uniprot, text = 'Orthologs'):
  uniprot = str(uniprot)
  return make_link(text, uniprot_orthologs_url(uniprot), 'orthologs') 

def uniprot_ppi_url(id):
  if id:
    return '/ppi/%s/' % id
  else:
    return ''

def uniprot_ppi_link(uniprot_id, text = 'PPI'):
  return make_link(text, uniprot_ppi_url(uniprot_id), 'ppi') 

"""
  Functions to parse Pfam designations from pfam_color_file and create links
"""

def make_pfam_list(family_name):
  handle = open_pfam_color_file(family_name)
  if handle:
    lines = handle.readlines()
    handle.close()
    return [line.split('=')[0] for line in lines]
  return []
    
def make_pfam_links(family_name, include_nonlinks = False):
  pfam_links = []
  handle = open_pfam_color_file(family_name)
  if handle:
    lines = handle.readlines()
    handle.close()
    pfams = [line.split('=')[0] for line in lines]
    tm_count = 0;
    for pfam in pfams:
      if pfam == 'TM':
        tm_count += 1;
    pfams = remove_duplicates(pfams)
    TM_found = SP_found = False
    for pfam in pfams:
      if not pfam in ['TM', 'SP']:
        pfam_links.append(pfam_link(pfam))
      else:
        TM_found = pfam == 'TM'
        SP_found = pfam == 'SP'
    # If included, force these non-links to the end of the list   
    if include_nonlinks: 
      if TM_found:
        if tm_count == 1:
          pfam_links.append('Transmembrane helix')
        else:
          pfam_links.append('Transmembrane helices [%s]' % tm_count)
      if SP_found:
        pfam_links.append('Signal peptide')
  return mark_safe(', '.join(pfam_links))

