#!/usr/bin/env python

from Bio import SwissProt
from Bio.SeqUtils import CheckSum
from optparse import OptionParser
from pfacts003.phylofacts.models import Sequence, SequenceHeader, \
    UniProtTaxonomy, UniProt, UniProtDatIndex, EC, UniProtEC, \
    GO_Term, GO_EvidencePriority, UniProtGO, PDB_Chain, PDB, \
    UniProtPDB_Chain, Pfam, UniProtPfam, UniProtHostOrganism, \
    UniProtFeature, UniProtGeneID, PostTranslationalModificationType, \
    FeatureKey, NonExperimentalQualifier, UniProtKeyword, Keyword, \
    UniProtOrganelle, Organelle, UniProtLiterature
import re



def handle_missing_taxonomy(record):

    """
    Handle rare case where UniProt taxonomy entry didn't exist

    Although one would expect all UniProt entries to have the
    prerequisite UniProt taxonomy entries before they are published,
    this does not always happen. We have seen several instances where
    the UniProt data has a reference to UniProt taxonomies that are not
    yet in the UniProt taxonomy database.

    For example, at the time of this writing, we see the following
    UniProt entry:

    ID   C8YT98_TAUSI            Unreviewed;        76 AA.
    AC   C8YT98;
    [snip]
    OS   Taumantis sigiana (Lime green praying mantis) (Taumantis ehr...
    OC   Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota;
    OC   Neoptera; Orthopteroidea; Dictyoptera; Mantodea; Mantidae;
    OC   Miomantinae; Miomantini; Taumantis.
    OX   NCBI_TaxID=765348;
    [snip]

    And, this taxon does not exist via the UniProt taxonomy page (where
    we retrieve our download data):
    http://www.uniprot.org/taxonomy/?query=765348&sort=score

    Since we have that NCBI TaxID and the taxon information, this data
    is inserted into the database so we can continue loading the UniProt
    data.

    It is believed that this entry will eventually get over-written by
    the correct entry whenever the UniProt taxonomy database is updated
    and we download and install it within our database.

    In case errors occur in the future, related to this problem, the
    other_names field also has the string "(Adding missing taxonomy)"
    attached. This should give us the flexibility to load the data, flag
    the data, and continue.

    We also print a message to standard out alerting the user that this
    has occurred. Standard out, instead of standard error, was used as
    this is not an uncaught error - but a warning for the user.

    The OS field is given to us in several formats. Examples of a
    typical OS line follow:

    OS   Escherichia coli.
    OS   Homo sapiens (Human).
    OS   Solanum melongena (Eggplant) (Aubergine).
    OS   Rous sarcoma virus (strain Schmidt-Ruppin A) (RSV-SRA) (Avian
    OS   leukosis virus-RSA).

    This is parsed into a single line by BioPython. As this is a
    temporary record, this entire string is inserted into the
    scientific_name to avoid any confusion with non-unique scientific
    names and incorrect parsing.

    After the taxonomy object is created, it is returned to the calling
    routine.
    """

    print "Warning: Taxonomy %s doesn't exist. Adding %s" % (
        record.taxonomy_id[0],
        record.organism)

    taxon = UniProtTaxonomy.objects.create(
        id=int(record.taxonomy_id[0]),
        scientific_name=record.organism,
        other_names="(Added missing taxonomy)"
    )
    return taxon

def main():
    parser = OptionParser(usage='%prog [StartPos_EndPos]')
    parser.add_option('--update_features', action='store_true',
      dest='update_features', default=True,
      help="Insert new entries into uniprot_feature and delete old ones.")
    parser.add_option('--no_update_features', action='store_false',
      dest='update_features', default=True,
      help="Only insert entries into uniprot_feature for new uniprot records.")
    (options, args) = parser.parse_args()
    sharding = False
    if len(args) >= 1:
        sharding = True
        shard_spec = args[0]
        positions = shard_spec.split('_')
        if positions < 2:
            parser.error(
              "Must specify shard as starting and ending file positions " \
              + "separated by an underscore")
        try:
            start_pos = int(positions[0])
        except ValueError:
            parser.error(
              "Must specify shard as starting and ending file positions " \
              + "separated by an underscore")
        try:
            end_pos = int(positions[1])
        except ValueError:
            parser.error(
              "Must specify shard as starting and ending file positions " \
              + "separated by an underscore")
    num_records = 0
    description_re = re.compile('(RecName: |AltName: |SubName: '
                + '|Full=|Short=|EC=|Allergen=|Biotech=|CD_antigen=|INN=|;)'
                + '|Includes: |Contains: |Flags: ')

    # Prepare regular expressions and object maps for parsing feature tables
    ptm_types = PostTranslationalModificationType.objects.all()
    ptm_re = re.compile('(%s)' % '|'.join([ptm_type.modification
                                            for ptm_type in ptm_types]))
    ptm_type_object_of_modification = {}
    for ptm_type in ptm_types:
        ptm_type_object_of_modification[ptm_type.modification] = ptm_type
    # Example:
    # RP   PHOSPHORYLATION [LARGE SCALE ANALYSIS] AT SER-267
    # Here PHOSPHORYLATION is a post-translational modification type, which
    # occurs at position 267.  Later there is a corresponding line in the
    # feature table:
    # FT   MOD_RES     267    267       Phosphoserine.
    ptm_pos_re = re.compile(' AT [A-Z][A-Z][A-Z]-([0-9]*)')
    feature_keys = FeatureKey.objects.all()
    feature_key_object_of_key_name = {}
    for feature_key in feature_keys:
        feature_key_object_of_key_name[feature_key.key_name] = feature_key
    nonexperimental_qualifiers = NonExperimentalQualifier.objects.all()
    nonexperimental_re = re.compile('(%s)' % '|'.join([qualifier.description
                                    for qualifier in nonexperimental_qualifiers]))
    nonexperimental_qualifier_object_of_description = {}
    for qualifier in nonexperimental_qualifiers:
        nonexperimental_qualifier_object_of_description[qualifier.description] \
            = qualifier
    dbSNPrs_re = re.compile('dbSNP:rs([0-9]*)')
    large_scale_re = re.compile('LARGE SCALE')

    go_evidence_objects = GO_EvidencePriority.objects.all()
    go_evidence_object_of_go_evidence_code = {}
    for go_evidence_object in go_evidence_objects:
        go_evidence_object_of_go_evidence_code[go_evidence_object.evidence]\
            = go_evidence_object

    f = open("/clusterfs/ohana/external/UniProt/to_import/uniprot.dat")
    pos_f = open("/clusterfs/ohana/external/UniProt/to_import/uniprot.dat")
    if sharding:
        f.seek(start_pos)
    current_pos = f.tell()
    pos_f.seek(current_pos)
    for record in SwissProt.parse(f):
        try:
            taxon = UniProtTaxonomy.objects.get(id__exact = record.taxonomy_id[0])
        except UniProtTaxonomy.DoesNotExist:
            taxon = handle_missing_taxonomy(record)
        seguid = CheckSum.seguid(record.sequence)

        # Parse the description
        description_tokens = description_re.split(record.description)
        full_recommended_name = ''
        # Look for the first recommended name category
        # (before any Includes or Contains sections)
        for i in xrange(len(description_tokens)):
            if description_tokens[i] == 'RecName: ':
                break
        # Now look for the full name
        for j in xrange(i, len(description_tokens)):
            if description_tokens[j] == 'Full=':
                break
        # The full recommended name is the next token
        if j < len(description_tokens) - 1:
            full_recommended_name = description_tokens[j+1]
        else:
            # Try looking for SubName instead, maybe this is a fragment
            # Look for the first subname category
            # (before any Includes or Contains sections)
            for i in xrange(len(description_tokens)):
                if description_tokens[i] == 'SubName: ':
                    break
            # Now look for the full name
            for j in xrange(i, len(description_tokens)):
                if description_tokens[j] == 'Full=':
                    break
            # The full subname is the next token
            if j < len(description_tokens) - 1:
                full_recommended_name = description_tokens[j+1]
            else:
                print "Full recommended name not found for %s" % record.entry_name
                print record.description
        # Look for all the EC numbers
        ecs = set()
        for i in xrange(len(description_tokens)):
            if description_tokens[i] == 'EC=':
                ecs.add(description_tokens[i+1])
        # Look for precursor or fragment flags
        is_fragment = False
        is_precursor = False
        for i in xrange(len(description_tokens)):
            if description_tokens[i] == 'Flag: ':
                if description_tokens[i+1][0:8] == 'Fragment':
                    is_fragment = True
                elif description_tokens[i+1] == 'Precursor':
                    is_precursor = True

        # Every UniProt accession is present in the uniprot_dat_index table.
        # Each of them points to a record in the uniprot table.
        # On the other hand, a record in the uniprot table has only one
        # accession, the one that was the primary accession the last time we did
        # this update.
        # The primary accession may have changed since we last updated (it may now
        # be a secondary accession).
        # The identifier may also have changed.  E.g., if the record was previously
        # in TrEMBL and is now in SwissProt, then its identifier may have changed
        # from one like Q197F8_IIV3 to one like 002R_IIV3 (i.e., the first part is
        # no longer the accession, but a gene name or something more informative).
        # So, we can't necessarily tell which was the existing record in the uniprot
        # table corresponding to the record we are now parsing by looking at either
        # its identifier or its accession.
        # Instead, we find the entries in the uniprot_dat_index table for each of
        # the accessions.
        # If one corresponding to the primary accession is already present, we take
        # the corresponding record in the UniProt table to be the canonical entry
        # corresponding to this UniProt record.
        # If none corresponding to the primary accession is present but entries in
        # the uniprot_dat_index for other accessions are present, we pick one of
        # these and make the corresponding record in the uniprot table the canonical
        # entry.
        # If no entries in the uniprot_dat_index table corresponding to any of these
        # accessions are present, we create a new record in the uniprot table and
        # make it the canonical entry.
        # If the entries in the uniprot_dat_index table corresponding to these
        # accessions point to multiple different records in the uniprot table, we
        # will delete the other ones at the end of this loop iteration.  But first
        # we will update the sequence_header records that point to those entries to
        # point instead to the canonical entry.
        new_uniprot = False
        uniprot_ids_to_delete = set()
        uniprot_dat_index_of_accession = {}
        uniprot_dat_indices = UniProtDatIndex.objects.filter(
                                uniprot_accession__in = record.accessions)
        uniprot_of_uniprot_id = {}
        for uniprot_dat_index in uniprot_dat_indices:
            uniprot_dat_index_of_accession[uniprot_dat_index.uniprot_accession] \
                = uniprot_dat_index
            uniprot_of_uniprot_id[uniprot_dat_index.uniprot.id] \
                = uniprot_dat_index.uniprot
        if len(uniprot_dat_index_of_accession.keys()) > 0:
            if record.accessions[0] in uniprot_dat_index_of_accession:
                uniprot = uniprot_dat_index_of_accession[record.accessions[0]].uniprot
            else:
                an_accession = uniprot_dat_index_of_accession.keys()[0]
                uniprot = uniprot_dat_index_of_accession[an_accession].uniprot
            for accession in uniprot_dat_index_of_accession:
                uniprot_dat_index_of_accession[accession].uniprot = uniprot
                uniprot_dat_index_of_accession[accession].file_char = current_pos
                uniprot_dat_index_of_accession[accession].save()
            missing_accessions \
              = set(record.accessions) - set(uniprot_dat_index_of_accession.keys())
            for accession in missing_accessions:
                uniprot_dat_index_of_accession[accession] \
                    = UniProtDatIndex.objects.create(file_char = current_pos,
                                              uniprot_accession = accession,
                                              uniprot = uniprot)
            uniprot_ids_to_delete = set(uniprot_of_uniprot_id.keys())
            uniprot_ids_to_delete.remove(uniprot.id)
            # Find sequence_headers pointing to the obsolete uniprot records, and
            # point them at the canonical one instead
            sequence_headers = SequenceHeader.objects.filter(uniprot__id__in =
                                                            uniprot_ids_to_delete)
            for sequence_header in sequence_headers:
                sequence_header.uniprot = uniprot
                sequence_header.save()
        else:
            new_uniprot = True
            uniprot = UniProt.objects.create(uniprot_identifier = record.entry_name,
                          accession = record.accessions[0],
                          taxon = taxon,
                          de = full_recommended_name,
                          seguid = seguid,
                          in_swissprot_f = (record.data_class == 'Reviewed'))
            for accession in record.accessions:
                uniprot_dat_index_of_accession[accession] \
                    = UniProtDatIndex.objects.create(file_char = current_pos,
                                              uniprot_accession = accession,
                                              uniprot = uniprot)
            # Look for orphaned sequence headers that can be assigned to this
            # uniprot
            sequences = Sequence.objects.filter(seguid__exact = seguid)
            sequence_ids = [sequence.id for sequence in sequences]
            if sequences:
                sequence_headers = SequenceHeader.objects.filter(
                                          sequence__id__in = sequence_ids,
                                          uniprot__isnull = True,
                                          taxon__id__exact = taxon.id)
                for sequence_header in sequence_headers:
                    sequence_header.uniprot = uniprot

        uniprot.uniprot_identifier = record.entry_name
        uniprot.accession = record.accessions[0]
        uniprot.uniprot_taxon = taxon
        uniprot.de = full_recommended_name
        uniprot.seguid = seguid
        uniprot.in_swissprot_f = (record.data_class == 'Reviewed')
        uniprot.description = record.description
        if is_fragment:
            uniprot.is_fragment = True
        if is_precursor:
            uniprot.is_precursor = True
        uniprot.save()

        # Update the EC associations
        uniprot_ec_objects = UniProtEC.objects.filter(uniprot__exact =
                                                              uniprot)
        uniprot_ec_object_of_ec_id = {}
        for uniprot_ec_object in uniprot_ec_objects:
            uniprot_ec_object_of_ec_id[uniprot_ec_object.ec.id] = uniprot_ec_object
        db_ec_ids = set(uniprot_ec_object_of_ec_id.keys())
        ec_object_of_ec_id = {}
        for ec in ecs:
            class_number_str, subclass_number_str, subsubclass_number_str, \
                enzyme_number_str = ec.split('.')
            is_preliminary = False

            # If EC number is similar to '-.-.-.-', then we should report
            # an error and continue. This is obviously dirty data and should
            # be reported to UniProt
            try:
                class_number = int(class_number_str)
            except ValueError:
                print "Warning: %s has invalid EC Number: '%s'" % (
                    uniprot.uniprot_identifier, ec)
                continue

            if subclass_number_str == '-':
                ec_objects = EC.objects.filter(class_number__exact = class_number,
                                          subclass_number__isnull = True,
                                          subsubclass_number__isnull = True,
                                          enzyme_number__isnull = True)
            else:
                subclass_number = int(subclass_number_str)
                if subsubclass_number_str == '-':
                    ec_objects = EC.objects.filter(class_number__exact = class_number,
                                              subclass_number__exact = subclass_number,
                                              subsubclass_number__isnull = True,
                                              enzyme_number__isnull = True)
                else:
                    subsubclass_number = int(subsubclass_number_str)
                    enzyme_number_str = enzyme_number_str.strip().rstrip(';')
                    if enzyme_number_str == '-':
                        ec_objects = EC.objects.filter(class_number__exact = class_number,
                                          subclass_number__exact = subclass_number,
                                          subsubclass_number__exact = subsubclass_number,
                                          enzyme_number__isnull = True)
                    else:
                        try:
                            enzyme_number = int(enzyme_number_str)
                        except ValueError:
                            print "Preliminary EC %s in %s" % (ec, record.entry_name)
                            print record.description
                            is_preliminary = True
                            enzyme_number = int(enzyme_number_str[1:])
                        ec_objects = EC.objects.filter(class_number__exact = class_number,
                                          subclass_number__exact = subclass_number,
                                          subsubclass_number__exact = subsubclass_number,
                                          enzyme_number__exact = enzyme_number,
                                          is_preliminary_f = is_preliminary)
            if ec_objects:
                ec_object = ec_objects[0]
            else:
                ec_object = EC.objects.create(class_number = class_number,
                                              subclass_number = subclass_number,
                                              subsubclass_number = subsubclass_number,
                                              enzyme_number = enzyme_number,
                                              is_preliminary_f = is_preliminary)
            ec_object_of_ec_id[ec_object.id] = ec_object
        uniprot_dat_ec_ids = set(ec_object_of_ec_id.keys())
        for ec_id in db_ec_ids - uniprot_dat_ec_ids:
            uniprot_ec_object_of_ec_id[ec_id].delete()
        for ec_id in uniprot_dat_ec_ids - db_ec_ids:
            UniProtEC.objects.create(uniprot = uniprot,
                                      ec = ec_object_of_ec_id[ec_id])

        # Update the keyword associations
        uniprot_keyword_objects = UniProtKeyword.objects.filter(uniprot__exact =
                                                              uniprot)
        uniprot_keyword_object_of_keyword_accession = {}
        for uniprot_keyword_object in uniprot_keyword_objects:
            uniprot_keyword_object_of_keyword_accession[ \
                  uniprot_keyword_object.keyword.accession] = uniprot_keyword_object
        db_keyword_accessions \
            = set(uniprot_keyword_object_of_keyword_accession.keys())
        keyword_object_of_keyword_accession = {}
        for keyword in record.keywords:
            keyword_objects = Keyword.objects.filter(identifier__exact = keyword)
            if keyword_objects:
                keyword_object = keyword_objects[0]
                keyword_object_of_keyword_accession[keyword_object.accession] \
                    = keyword_object
            else:
                print "Unrecognized keyword %s while parsing %s." % (keyword,
                    record.entry_name),
                print "The keyword table may be out of date."
        uniprot_dat_keyword_accessions \
            = set(keyword_object_of_keyword_accession.keys())
        for accession in db_keyword_accessions - uniprot_dat_keyword_accessions:
            uniprot_keyword_object_of_keyword_accession[accession].delete()
        for accession in uniprot_dat_keyword_accessions - db_keyword_accessions:
            UniProtKeyword.objects.create(uniprot = uniprot,
                    keyword = keyword_object_of_keyword_accession[accession])

        # Update the organelle associations
        uniprot_organelle_objects = UniProtOrganelle.objects.filter(
                                                        uniprot__exact = uniprot)
        uniprot_organelle_object_of_organelle_id = {}
        for uniprot_organelle_object in uniprot_organelle_objects:
            uniprot_organelle_object_of_organelle_id[ \
                uniprot_organelle_object.organelle.id] = uniprot_organelle_object
        db_organelle_ids = set(uniprot_organelle_object_of_organelle_id.keys())
        organelle_object_of_organelle_id = {}
        for organelle in record.organelle.rstrip('.').split(','):
            if len(organelle) == 0:
                continue
            fields = organelle.split('; ')
            if len(fields) > 1:
                # This had better be a plastid
                if fields[0] == 'Plastid':
                    organelle_objects = Organelle.objects.filter(
                                            description__exact = fields[0],
                                            plastid_type__exact = fields[1])
                else:
                    print "Unrecognized organelle %s in %s" % (organelle,
                                                                record.entry_name)
            else:
                organelle_objects = Organelle.objects.filter(
                                        description__exact = fields[0])
            if organelle_objects:
                organelle_object = organelle_objects[0]
            else:
                field = fields[0].strip()
                # This had better be a plasmid
                if len(field) >= 9 and field[0:7] == 'Plasmid':
                    organelle_object = Organelle.objects.create(
                                                  description = field,
                                                  plasmid_name = field[8:])
                elif len(field) >= 13 and field[4:11] == 'Plasmid':
                    organelle_object = Organelle.objects.create(
                                                  description = field[4:],
                                                  plasmid_name = field[12:])
                else:
                    print "Unrecognized organelle %s in %s" % (organelle,
                                                                record.entry_name)
                    continue
            organelle_object_of_organelle_id[organelle_object.id] = organelle_object
        uniprot_dat_organelle_ids = set(organelle_object_of_organelle_id.keys())
        for organelle_id in db_organelle_ids - uniprot_dat_organelle_ids:
            uniprot_organelle_object_of_organelle_id[organelle_id].delete()
        for organelle_id in uniprot_dat_organelle_ids - db_organelle_ids:
            UniProtOrganelle.objects.create(uniprot = uniprot,
                        organelle = organelle_object_of_organelle_id[organelle_id])

        # Update the host organism associations
        uniprot_host_objects = UniProtHostOrganism.objects.filter(uniprot__exact =
                                                              uniprot)
        uniprot_host_object_of_host_id = {}
        for uniprot_host_object in uniprot_host_objects:
            uniprot_host_object_of_host_id[uniprot_host_object.host_organism.id] \
                = uniprot_host_object
        db_host_ids = set(uniprot_host_object_of_host_id.keys())
        host_object_of_host_id = {}
        host_ids = [int(host_spec.split(';')[0])
                    for host_spec in record.host_organism]
        host_objects = UniProtTaxonomy.objects.filter(id__in = host_ids)
        for host_object in host_objects:
            host_object_of_host_id[host_object.id] = host_object
        uniprot_dat_host_ids = set(host_object_of_host_id.keys())
        for host_id in set(host_ids) - uniprot_dat_host_ids:
            print "Unknown host taxonomy id %d when parsing %s;" \
                % (host_id, record.entry_name),
            print "UniProtTaxonomy table may be out of date."
        for host_id in db_host_ids - uniprot_dat_host_ids:
            uniprot_host_object_of_host_id[host_id].delete()
        for host_id in uniprot_dat_host_ids - db_host_ids:
            UniProtHostOrganism.objects.create(uniprot = uniprot,
                                    host_organism = host_object_of_host_id[host_id])

        # Update the literature references
        uniprot_literature_objects = UniProtLiterature.objects.filter(
                                                          uniprot__exact = uniprot)
        uniprot_literature_object_of_title = {}
        for uniprot_literature_object in uniprot_literature_objects:
            uniprot_literature_object_of_title[uniprot_literature_object.title] \
                = uniprot_literature_object
        db_titles = set(uniprot_literature_object_of_title.keys())
        uniprot_dat_ref_of_title = {}
        for ref in record.references:
            titles = ref.title.split(';')
            for title in titles:
                uniprot_dat_ref_of_title[title.strip('"')] = ref
        uniprot_dat_titles = set(uniprot_dat_ref_of_title.keys())
        for title in db_titles - uniprot_dat_titles:
            uniprot_literature_object_of_title[title].delete()
        for title in uniprot_dat_titles:
            if title in db_titles:
                uniprot_literature_object = uniprot_literature_object_of_title[title]
            else:
                uniprot_literature_object \
                  = UniProtLiterature.objects.create(uniprot = uniprot, title = title)
            ref = uniprot_dat_ref_of_title[title]
            uniprot_literature_object.authors = ref.authors
            positional_info = ' '.join(ref.positions)
            m = large_scale_re.search(positional_info)
            if m:
                uniprot_literature_object.is_large_scale_f = True
            for db_name, db_reference in ref.references:
                if db_name == 'MEDLINE':
                    uniprot_literature_object.medline_ui = db_reference
                elif db_name == 'PubMed':
                    uniprot_literature_object.pmid = db_reference
                elif db_name == 'DOI':
                    uniprot_literature_object.doi = db_reference
                elif db_name == 'AGRICOLA':
                    uniprot_literature_object.agricola = db_reference
            uniprot_literature_object.save()

        # Find cross references to other databases
        geneids = set()
        go_evidence_of_go_accession= {}
        pfam_accessions = set()
        extent_of_pdb_chain_ids = {}
        for reference in record.cross_references:
            if reference[0] == 'GeneID':
                geneids.add(reference[1])
            elif reference[0] == 'GO':
                if len(reference) >= 4:
                    go_evidence_of_go_accession[reference[1]] \
                        = reference[3].split(':')[0]
            elif reference[0] == 'Pfam':
                pfam_accessions.add(reference[1])
            elif reference[0] == 'PDB':
                pdb_id = reference[1].lower()
                chain_ids = reference[4].split('=')[0].split('/')
                pdb_from_residue = None
                pdb_to_residue = None
                fields = reference[4].split('=')
                if len(fields) > 1:
                    try:
                        pdb_from_residue, pdb_to_residue \
                            = [int(x) for x in fields[1].split('-')]
                    except IndexError:
                        pdb_from_residue = None
                        pdb_to_residue = None
                    except ValueError:
                        pdb_from_residue = None
                        pdb_to_residue = None
                for chain_id in chain_ids:
                    pdb_chain_id = pdb_id + chain_id
                    extent_of_pdb_chain_ids[pdb_chain_id] \
                        = (pdb_from_residue, pdb_to_residue)

        # Update the GeneID associations
        uniprot_geneid_objects = UniProtGeneID.objects.filter(uniprot__exact =
                                                              uniprot)
        uniprot_geneid_object_of_geneid = {}
        for uniprot_geneid_object in uniprot_geneid_objects:
            uniprot_geneid_object_of_geneid[uniprot_geneid_object.geneid] \
                = uniprot_geneid_object
        db_geneids = set(uniprot_geneid_object_of_geneid.keys())
        for geneid in db_geneids - geneids:
            uniprot_geneid_object_of_geneid[geneid].delete()
        for geneid in geneids - db_geneids:
            UniProtGeneID.objects.create(uniprot = uniprot, geneid = geneid)

        # Update the GO associations
        uniprot_go_objects = UniProtGO.objects.filter(uniprot__exact = uniprot)
        uniprot_go_object_of_go_term_accession = {}
        for uniprot_go_object in uniprot_go_objects:
            uniprot_go_object_of_go_term_accession[uniprot_go_object.go_term.acc] \
                = uniprot_go_object
        db_go_term_accessions = set(uniprot_go_object_of_go_term_accession.keys())
        go_term_objects = GO_Term.objects.filter(acc__in =
                                            go_evidence_of_go_accession.keys())
        go_term_object_of_go_term_accession = {}
        for go_term_object in go_term_objects:
            go_term_object_of_go_term_accession[go_term_object.acc] \
                = go_term_object
        uniprot_dat_go_term_accessions \
            = set(go_term_object_of_go_term_accession.keys())
        for go_accession in set(go_evidence_of_go_accession.keys()) - \
                              uniprot_dat_go_term_accessions:
            print "Unrecognized GO accession %s while parsing %s" \
                % (go_accession, record.entry_name),
            print "GO term table may be out of date."
        for go_evidence_code in set(go_evidence_of_go_accession.values()) - \
                          set(go_evidence_object_of_go_evidence_code.keys()):
            print "Unrecognized GO evidence code %s while parsing %s" \
                % (go_evidence_code, record.entry_name),
            print "GO evidence_priority table may be out of date."
        for go_term_accession in db_go_term_accessions \
                                  - uniprot_dat_go_term_accessions:
            uniprot_go_object_of_go_term_accession[go_term_accession].delete()
        for go_term_accession in uniprot_dat_go_term_accessions:
            go_evidence_code = go_evidence_of_go_accession[go_term_accession]
            if go_evidence_code in go_evidence_object_of_go_evidence_code:
                go_evidence_object = \
                    go_evidence_object_of_go_evidence_code[go_evidence_code]
                if go_term_accession in db_go_term_accessions:
                    uniprot_go \
                        = uniprot_go_object_of_go_term_accession[go_term_accession]
                    if uniprot_go.go_evidence.evidence != go_evidence_code:
                        uniprot_go.go_evidence = go_evidence_object
                        uniprot_go.save()
                else:
                    UniProtGO.objects.create(go_term =
                                go_term_object_of_go_term_accession[go_term_accession],
                                go_evidence = go_evidence_object,
                                uniprot = uniprot)

        # Update the Pfam associations
        uniprot_pfam_objects = UniProtPfam.objects.filter(uniprot__exact = uniprot)
        uniprot_pfam_object_of_pfam_accession = {}
        for uniprot_pfam_object in uniprot_pfam_objects:
            uniprot_pfam_object_of_pfam_accession[
                          uniprot_pfam_object.pfam.accession] = uniprot_pfam_object
        db_pfam_accessions = set(uniprot_pfam_object_of_pfam_accession.keys())
        pfam_object_of_pfam_accession = {}
        for pfam_accession in pfam_accessions:
            pfam_objects = Pfam.objects.filter(
                                        accession__exact=pfam_accession).order_by(
                                        'overall_pfam_version').reverse()
            if pfam_objects:
                pfam_object = pfam_objects[0]
                pfam_object_of_pfam_accession[pfam_object.accession] = pfam_object
            else:
                print "Unknown Pfam accession %s encountered when parsing %s" \
                    % (pfam_accession, record.entry_name),
                print "Pfam table may be out of date"
        uniprot_dat_pfam_accessions = set(pfam_object_of_pfam_accession.keys())
        for pfam_accession in db_pfam_accessions - uniprot_dat_pfam_accessions:
            uniprot_pfam_object_of_pfam_accession[pfam_accession].delete()
        for pfam_accession in uniprot_dat_pfam_accessions - db_pfam_accessions:
            UniProtPfam.objects.create(uniprot = uniprot,
                              pfam = pfam_object_of_pfam_accession[pfam_accession])

        # Update the PDB associations
        uniprot_pdb_chain_objects = UniProtPDB_Chain.objects.filter(uniprot__exact =
                                                                    uniprot)
        uniprot_pdb_chain_object_of_pdb_chain_id = {}
        for uniprot_pdb_chain_object in uniprot_pdb_chain_objects:
            pdb_chain_id = uniprot_pdb_chain_object.pdb_chain.pdb.id + \
                            uniprot_pdb_chain_object.pdb_chain.chain_id
            uniprot_pdb_chain_object_of_pdb_chain_id[pdb_chain_id] \
                = uniprot_pdb_chain_object
        db_pdb_chain_ids = set(uniprot_pdb_chain_object_of_pdb_chain_id.keys())
        pdb_chain_object_of_pdb_chain_id = {}
        for pdb_chain_id in extent_of_pdb_chain_ids.keys():
            pdb_chain_objects = PDB_Chain.objects.filter(
                                      pdb__id__exact = pdb_chain_id[0:4],
                                      chain_id__exact = pdb_chain_id[4:])
            if pdb_chain_objects:
                pdb_chain_object = pdb_chain_objects[0]
                pdb_chain_object_of_pdb_chain_id[pdb_chain_id] = pdb_chain_object
            else:
                print "Unknown PDB chain %s encountered when parsing %s" \
                    % (pdb_chain_id, record.entry_name),
                print "The PDB_Chain table may be out of date."
        uniprot_dat_pdb_chain_ids = set(pdb_chain_object_of_pdb_chain_id.keys())
        for pdb_chain_id in db_pdb_chain_ids - uniprot_dat_pdb_chain_ids:
            uniprot_pdb_chain_object_of_pdb_chain_id[pdb_chain_id].delete()
        for pdb_chain_id in uniprot_dat_pdb_chain_ids:
            pdb_from_residue, pdb_to_residue = extent_of_pdb_chain_ids[pdb_chain_id]
            if pdb_from_residue:
                if pdb_chain_id in db_pdb_chain_ids:
                    uniprot_pdb_chain \
                        = uniprot_pdb_chain_object_of_pdb_chain_id[pdb_chain_id]
                    uniprot_pdb_chain.from_residue = pdb_from_residue
                    uniprot_pdb_chain.to_residue = pdb_to_residue
                    uniprot_pdb_chain.save()
                else:
                    UniProtPDB_Chain.objects.create(uniprot = uniprot,
                              pdb_chain = pdb_chain_object_of_pdb_chain_id[pdb_chain_id],
                              from_residue = pdb_from_residue,
                              to_residue = pdb_to_residue)
            else:
                UniProtPDB_Chain.objects.create(uniprot = uniprot,
                            pdb_chain = pdb_chain_object_of_pdb_chain_id[pdb_chain_id])

        # Update the feature table (position-specific information)
        if new_uniprot or options.update_features:
            # Unfortunately, there is no part of a feature table entry that is
            # guaranteed to persist from one release of UniProt of the next (except
            # the FTid, but that's not always present).  Therefore there's no way to
            # easily determine that a feature in the uniprot.dat line is the same or
            # nearly the same as one that's already in the database.  So, we update
            # the features by simply inserting all the entries anew and deleting the
            # old ones (without trying to check if they were the same or modify the
            # old ones).  So, we get the old feature entries first so we can delete
            # them later.
            uniprot_feature_objects = UniProtFeature.objects.filter(
                                                          uniprot__exact = uniprot)
            # Instantiate the queryset by turning it into a list.  Otherwise, it
            # won't be instantiated until we get to the bottom and loop over these to
            # delete them, at which point it would delete *all* of them (including
            # the ones we just created).
            uniprot_feature_object_list = list(uniprot_feature_objects)

            # The fields from_residue_is_uncertain, to_residue_is_uncertain,
            # extends_n_terminally, and extends_c_terminally derive from the FT line
            # in the uniprot.dat line, according to the UniProt KnowledgeBase user
            # manual:

            # When a feature is known to extend beyond the position that is given in
            # the feature table, the endpoint specification will be preceded by '<'
            # for features which continue to the left end (N-terminal direction) or
            # by '>' for features which continue to the right end (C- terminal
            # direction); Unknown endpoints are denoted by '?'. Uncertain endpoints
            # are denoted by a '?' before the position, e.g. '?42'.

            for feature in record.features:
                key_name, from_residue_spec, to_residue_spec, description, \
                    ftid = feature
                if key_name in feature_key_object_of_key_name:
                    feature_key = feature_key_object_of_key_name[key_name]
                    # Check for nonexperimental qualifier
                    nonexperimental_qualifier = None
                    m = nonexperimental_re.search(description)
                    if m:
                        nonexperimental_qualifier \
                            = nonexperimental_qualifier_object_of_description[m.group(0)]
                    created_object = False
                    if key_name == 'VARIANT':
                        # Check for dbSNP:rsaccession_number
                        dbSNP_rs_accession = None
                        m = dbSNPrs_re.search(description)
                        if m:
                            dbSNP_rs_accession = int(m.group(1))
                            if nonexperimental_qualifier:
                                if ftid == '':
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          dbsnp_rs_number = dbSNP_rs_accession,
                                          nonexperimental_qualifier = nonexperimental_qualifier)
                                else:
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          feature_identifier = ftid,
                                          dbsnp_rs_number = dbSNP_rs_accession,
                                          nonexperimental_qualifier = nonexperimental_qualifier)
                            else:
                                if ftid == '':
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          dbsnp_rs_number = dbSNP_rs_accession)
                                else:
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          feature_identifier = ftid,
                                          dbsnp_rs_number = dbSNP_rs_accession)
                            created_object = True
                    elif key_name == 'MOD_RES' and from_residue_spec == to_residue_spec:
                        # Look for the post-translational modification type
                        found_ptm = False
                        for reference in record.references:
                            if found_ptm:
                                break
                            positional_info = ' '.join(reference.positions)
                            ptm_tokens = ptm_re.split(positional_info)
                            if len(ptm_tokens) > 1:
                                for i in range((len(ptm_tokens)-1) / 2):
                                    modification = ptm_tokens[2*i+1]
                                    string_with_position = ptm_tokens[2*i+2]
                                    match_position = ptm_pos_re.search(string_with_position)
                                    if match_position:
                                        position = int(match_position.group(1))
                                        if position == from_residue:
                                            # Success!
                                            found_ptm = True
                                            ptm_type = ptm_type_object_of_modification[modification]
                                            break
                        if found_ptm:
                            if nonexperimental_qualifier:
                                if ftid == '':
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          posttranslational_modification_type = ptm_type,
                                          nonexperimental_qualifier = nonexperimental_qualifier)
                                else:
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          feature_identifier = ftid,
                                          posttranslational_modification_type = ptm_type,
                                          nonexperimental_qualifier = nonexperimental_qualifier)
                            else:
                                if ftid == '':
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          posttranslational_modification_type = ptm_type)
                                else:
                                    feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                          feature_key = feature_key,
                                          description = description,
                                          feature_identifier = ftid,
                                          posttranslational_modification_type = ptm_type)
                            created_object = True
                    if not created_object:
                        if nonexperimental_qualifier:
                            if ftid == '':
                                feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                        feature_key = feature_key,
                                        description = description,
                                        nonexperimental_qualifier = nonexperimental_qualifier)
                            else:
                                feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                        feature_key = feature_key,
                                        description = description,
                                        feature_identifier = ftid,
                                        nonexperimental_qualifier = nonexperimental_qualifier)
                        else:
                            if ftid == '':
                                feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                        feature_key = feature_key,
                                        description = description)
                            else:
                                feature_obj = UniProtFeature.objects.create(uniprot = uniprot,
                                        feature_key = feature_key,
                                        description = description,
                                        feature_identifier = ftid)
                else:
                    print "Unrecognized feature key %s while parsing %s" \
                        % (key_name, record.entry_name)
                # Parse the from_residue and to_residue and update the appropriate
                # fields in the object.
                try:
                    from_residue = int(from_residue_spec)
                    feature_obj.from_residue = from_residue
                    feature_obj.save()
                except ValueError:
                    if from_residue_spec != '?':
                        for i in range(len(from_residue_spec)):
                            if from_residue_spec[i] == '<':
                                feature_obj.extends_n_terminally = True
                            elif from_residue_spec[i] == '>':
                                # We don't expect this to happen, but anyway...
                                feature_obj.extends_c_terminally = True
                            elif from_residue_spec[i] == '?':
                                feature_obj.from_residue_is_uncertain = True
                            else:
                                break
                        feature_obj.from_residue = int(from_residue_spec[i:])
                        feature_obj.save()
                try:
                    to_residue = int(to_residue_spec)
                    feature_obj.to_residue = to_residue
                    feature_obj.save()
                except ValueError:
                    if to_residue_spec != '?':
                        for i in range(len(to_residue_spec)):
                            if to_residue_spec[i] == '<':
                                # We don't expect this to happen, but anyway...
                                feature_obj.extends_n_terminally = True
                            elif to_residue_spec[i] == '>':
                                feature_obj.extends_c_terminally = True
                            elif to_residue_spec[i] == '?':
                                feature_obj.to_residue_is_uncertain = True
                            else:
                                break
                        feature_obj.to_residue = int(to_residue_spec[i:])
                        feature_obj.save()
            # Delete the old feature entries.
            for uniprot_feature_object in uniprot_feature_object_list:
                uniprot_feature_object.delete()


        for id in uniprot_ids_to_delete:
            uniprot_of_uniprot_id[id].delete()
        # The SwissProt parser may have eaten many characters from the next record
        # by this point, buffering them away until we ask for the next record.  So
        # if we set current_pos from f.tell() now, we will get the wrong answer.
        # Instead, we will update current_pos by reading lines from pos_f, without
        # parsing them--we're just looking for the record separator.
        line = pos_f.readline()
        while len(line) >= 2 and line[0:2] != '//':
            line = pos_f.readline()
        current_pos = pos_f.tell()
        num_records += 1
        if sharding:
            if current_pos > end_pos:
                break

if __name__ == '__main__':
    main()
