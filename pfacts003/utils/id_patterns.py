import re

uniprot_accession_pat1 = '[OoPpQq][0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][0-9]'
uniprot_accession_pat2 \
  = '[A-NR-Za-nr-z][0-9][A-Za-z][A-Za-z0-9][A-Za-z0-9][0-9]'
uniprot_taxon_pat \
  = '[A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9]?[A-Za-z0-9]?'
swissprot_desc_pat = \
  '[A-Za-z0-9][A-Za-z0-9]?[A-Za-z0-9]?[A-Za-z0-9]?[A-Za-z0-9]?'
uniprot_accession_re1 = re.compile('^%s$' % uniprot_accession_pat1)
uniprot_accession_re2 = re.compile('^%s$' % uniprot_accession_pat2)
uniprot_identifier_re1 = re.compile('^%s_%s$' % (uniprot_accession_pat1,
                                                uniprot_taxon_pat))
uniprot_identifier_re2 = re.compile('^%s_%s$' % (uniprot_accession_pat2,
                                                uniprot_taxon_pat))
uniprot_identifier_re3 = re.compile('^%s_%s$' % (swissprot_desc_pat,
                                                uniprot_taxon_pat))
gi_pat = '[0-9]+'
gi_re = re.compile(gi_pat)
taxon_id_pat = gi_pat
taxon_id_re = gi_re
leftid_pat = gi_pat
leftid_re = re.compile(leftid_pat)
aa_pat = '[AC-IK-NP-TVWY]*'
aa_with_whitespace_pat = '[AC-IK-NP-TVWY\n\t ]*'
fasta_pat = '(>[^\n]*\n)?(%s)\*?' % aa_with_whitespace_pat
aa_re = re.compile(aa_pat)
fasta_re = re.compile(fasta_pat)
scopid_pat = '[dg][0-9][0-9a-z][0-9a-z][0-9a-z_][0-9a-z_\.][1-9a-z_]'
scopid_re = re.compile(scopid_pat)
bpgid_pat = 'bpg[0-9][0-9][0-9][0-9][0-9][0-9][0-9]'
bpgid_re = re.compile(bpgid_pat)
uncharacterized_pat = '[Uu]ncharacterized [Pp]rotein'
uncharacterized_re = re.compile(uncharacterized_pat)
hypothetical_pat = '[Hh]ypothetical [Pp]rotein'
hypothetical_re = re.compile(hypothetical_pat)
predicted_pat = '[Pp]redicted [Pp]rotein'
predicted_re = re.compile(predicted_pat)
uncertain_re = re.compile(r'(uncharacterized|hypothetical|predicted) protein',
    re.I)
float_pat = '\d+(\.\d*)?([eE][-+]?\d+)?'
float_re = re.compile(float_pat)
# Float patterns that avoid the use of parentheses,
# since Django turns each parenthesized expression in urls.py into an argument
float_pat1 = '\d+'
float_pat2 = '\d+\.\d*'
float_pat3 = '\d+[eE][-+]?\d+'
float_pat4 = '\d+\.\d*[eE][-+]?\d+'
phog_pat = '[Pp][Hh][Oo][Gg]\d\d\d\d\d\d\d_\d\d\d\d\d'
phog_re = re.compile(phog_pat)
phog_tight_pat = '%sTC' % phog_pat
phog_tight_re = re.compile(phog_tight_pat)
phog_medium_pat = '%sTM' % phog_pat
phog_medium_re = re.compile(phog_medium_pat)
phog_loose_pat = '%sTD' % phog_pat
phog_loose_re = re.compile(phog_loose_pat)
phog_custom_pat = '%sT%s' % (phog_pat, float_pat)
phog_custom_re = re.compile(phog_custom_pat)
phog_custom_pat1 = '%sT%s' % (phog_pat, float_pat1)
phog_custom_pat2 = '%sT%s' % (phog_pat, float_pat2)
phog_custom_pat3 = '%sT%s' % (phog_pat, float_pat3)
phog_custom_pat4 = '%sT%s' % (phog_pat, float_pat4)
kegg_map_pat = 'map\d\d\d\d\d'
kegg_map_re = re.compile(kegg_map_pat)


def is_uniprot_identifier_format(alleged_identifier):
    """Verify if alleged uniprot identifier has a valid format

    UniProt identifiers are typically a pattern of two items combined with an
    underscore, where the first item is a gene accession and the second is the
    taxonomic mnemonic (e.g., RUBR1_PSEAE). The exact patterns for the
    identifier vary and are defined above. 

    This helper function takes an alleged_identifier, interets its format, and
    then returns either True or False if it matches or not (respectively). This
    method does not validate if the UniProt identifier actually exists within
    the UniProt database.
    """

    is_match = False
    
    if uniprot_identifier_re1.match(alleged_identifier) is not None or \
       uniprot_identifier_re2.match(alleged_identifier) is not None or \
       uniprot_identifier_re3.match(alleged_identifier) is not None:
        is_match = True

    return is_match

def is_uniprot_accession_format(alleged_accession):
    """Verify if alleged uniprot accession has a valid format

    UniProt accessions (not to be confused with UniProt identifiers) are a
    string of six characters with a specific set of rules for matching. See:
    http://www.uniprot.org/manual/accession_numbers for more details.

    This helper function takes an alleged_accession, interets its format, and
    then returns either True or False if it matches or not (respectively). This
    method does not validate if the UniProt accession actually exists within
    the UniProt database.
    """

    is_match = False
    
    if uniprot_accession_re1.match(alleged_accession) is not None or \
       uniprot_accession_re2.match(alleged_accession) is not None:
        is_match = True

    return is_match
