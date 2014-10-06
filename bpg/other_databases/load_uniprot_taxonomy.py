#!/usr/bin/env python

import sys

from django.db import transaction

from pfacts003.phylofacts import models


@transaction.commit_on_success
def update_uniprot_taxonomy(filename):
    """Create or Update Uniprot Taxonomic Data

    If the taxonomic data already exists, then the record is retrieved
    and updated to the correct value. If the data does not exist, then
    it is created.

    When this routine finishes, all of the taxonomic data downloaded
    from uniprot should be in the database.

    However, what is *not* completed is any updates to the left and
    right ids (which could now be completely wrong).
    """
    handle = open(filename, 'r')
    new_species = []
    for line in handle.readlines()[1:]:
        # Replace blank fields with None
        record = [x != '' and x or None for x in line.strip('\n').split('\t')]
        # If the "Reviewed" field is None, replace it with False;
        # if the "Reviewed" field was a non-blank string, replace it with True
        record[6] = record[6] is not None and True or False
        # If the "Parent" field was None, replace it with 0

        u, created = models.UniProtTaxonomy.objects.get_or_create(id=record[0]);
        u.mnemonic=record[1]
        u.scientific_name=record[2]
        u.common_name=record[3]
        u.synonym=record[4]
        u.other_names=record[5]
        u.reviewed=record[6]
        u.rank=record[7]

        if record[9] is not None:
            u.parent_id=record[9]

        u.save()

        if created:
            new_species.append(u)

    print "%d new species were added:" % len(new_species)
    for specie in new_species:
        print specie.id, specie

if __name__=="__main__":
    update_uniprot_taxonomy('/clusterfs/ohana/external/taxonomy/to_import/taxonomy-all.tab')
