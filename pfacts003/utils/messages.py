from django.utils.safestring import mark_safe

unrecognized_identifier = mark_safe("""Unrecognized identifier/accession &ndash;
<a href="/orthologs/faq#unrecognized_id">Why is this?</a>""")

def not_in_phog(sequence_id, family=False):
   if family:
       container = 'PhyloFacts family'
   else:
       container = 'PHOG'
   return mark_safe("""%s is not contained in any %s &ndash; 
       <a href="/phog/orthologs/faq#not_in_phog">Why is this?</a>""" % (sequence_id, container))

default_error_message = mark_safe("There was an error processing your request")

family_accession_format_message = mark_safe("""Your input must begin with 
the letters bpg, followed by exactly 7 digits.""")

no_uniprot_record_found = mark_safe("""No record was found for the information 
you entered.  Please enter valid UniProt Accessions or Identifiers only.""")

invalid_acc_or_ident_message = mark_safe("""Your input must either be a 
valid UniProt Accession or Identifier.""")

family_does_not_exist = mark_safe("""A PhyloFacts Family does not exist for
the family accession provided.""")

