'''
This knows about our DB, and it knows how to load new ontologies, and replace annotations
for an ontology
'''
from pfacts003.phylofacts.models import Annotation, Ontology, OntologyTerm, \
        OntologySubset, OntologySubsetMembership, OntologyDBXRef, Evidence, \
        UniProt

from django.db import transaction, connection
#TODO closure table

@transaction.commit_on_success
def load_ontology(name="", subsets=[], terms=[], subset_membership=[],
           connections=[], xrefs=[]):
    #step one - remove current version
    o = Ontology.objects.filter(name=name)
    if o.count() > 0:
        o[0].delete()
    #step two - insert the new one
    #two a - insert into ontology
    o = Ontology.objects.create(name=name)
    #two b - insert into terms
    for term in terms:
        o.terms.create(
            accession = term[0],
            name = term[1]
        )
    #two c - insert into subsets
    for subset in subsets:
        o.subsets.create(
            name = subset[0],
        )
    #two d - insert into subset_membership
    for membership in subset_membership:
        OntologySubsetMembership.objects.create(
            term=OntologyTerm.objects.get(accession=membership[0]),
            subset=OntologySubset.objects.get(name=membership[1])
        )
    #two e - insert into connections
        #TODO
    #two f - insert into xrefs
    for xref in xrefs:
        OntologyDBXRef.objects.create(
            from_term = OntologyTerm.objects.get(accession=xref[0]),
            to_accession = xref[1],
            type = "TODO"
        )


#@transaction.commit_on_success
def load_annotations(annotations, source):
    Annotation.objects.filter(source=source).delete()

    from guppy import hpy
    import sys
    hp = hpy()
    hp.setrelheap()
    # create local dict of uniprot accession to id and
    # the sets of available codes and term accessions

    codes = set([ e['code'] for e in Evidence.objects.values('code') ])
    terms = set([ o['accession'] for o in OntologyTerm.objects.values('accession') ])
    uniprot_dict = dict((row['accession'], row['id']) for row in UniProt.objects.values('id', 'accession'))

    counter = 0

    cursor = connection.cursor()
    for annotation in annotations:
        counter += 1

        if counter % 10000 == 0:
            transaction.commit()
            print hp.heap()
            print "\n\n+++\n\n"

        if annotation[0] in uniprot_dict and annotation[1] in terms  and annotation[2] in codes:
            insert_string = "INSERT INTO annotation (uniprot_id, ontology_accession, source, evidence_code, assigned_by, date_assigned) " +\
                        " VALUES (%s, '%s', '%s', '%s', '%s', '%s' ) " % (uniprot_dict[annotation[0]], annotation[1], source, annotation[2], annotation[3], annotation[4]) 
            cursor.execute(insert_string)

        transaction.commit()


def create_annotations(annotations, source, outfile):
    codes = set([ e['code'] for e in Evidence.objects.values('code') ])
    terms = set([ o['accession'] for o in OntologyTerm.objects.values('accession') ])
    uniprot_dict = dict((row['accession'], row['id']) for row in UniProt.objects.values('id', 'accession'))

    for annotation in annotations:
        if annotation[0] in uniprot_dict and annotation[1] in terms  and annotation[2] in codes:
            insert_string = "%s,%s,%s,%s,%s,%s\n" % (uniprot_dict[annotation[0]], 
                    annotation[1], source, annotation[2], annotation[3], annotation[4])
            outfile.write(insert_string)
