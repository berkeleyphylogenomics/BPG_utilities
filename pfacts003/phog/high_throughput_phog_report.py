#!/usr/bin/env python


"""
High Throughput Ortholog Report Generation

Ortholog data of different thresholds has been stored within the
database. This module retrieves that data and summarizes is so the
report can be available for download.

Currently, the input to this module is a set of genomes by both JGI
accession and UniProt accession identifiers, separated by a single
space. The first column is the original ID and the second is the UniProt
identifier. This report retrieves data based upon the second column. An
example input follows:

637017740,O24856_HELPY
637017742,O24857_HELPY
637017757,O24868_HELPY
637017758,O24869_HELPY
637017764,O24873_HELPY
....

If only a single set of identifiers are avilable (e.g., only UniProt
identifiers are available), this report can still be used in its current
format by 'padding' dummy identifiers in the first column. For example:

123456789,O24856_HELPY
123456789,O24857_HELPY
123456789,O24868_HELPY

The final report that is generated is a comma separated report with the
excel dialect with the following columns:

* Input Identifier * UniProt Identifier * Annotation Type * Annotation
* Experimental Evidence? * PHOGs * Orthologs

This resultant CSV file has exactly seven columns, with one row for each
distinct annotation of an individual sequence.  However, each individual
sequence may have many different annotations.  Therefore, a single
sequence may be represented by many rows in the table.  In the example
below, there are 3 rows corresponding to 3 annotations for the sequence
6371892, two rows corresponding to two annotations for the sequence
6371893, and 2 rows corresponding to 2 annotations for the sequence
6371894.

The advantage of this format is that less data is required to fit into
a single column in the case where a single sequence has many different
predicted GO functions, GO processes, EC numbers, etc.  Critically, this
format allows us to store the provenance of each annotation. This table
can be used to derive summary data for each gene.

For example:

6371892,O36572_HELPY,EC,2.3.1.13,N,"PHOG0103256,PHOG0372856","O26738_...
6371892,O36572_HELPY,KEGG,map027837,N,"PHOG0103256,PHOG0372856","O267...
6371892,O36752_HELPY,EC,2.3.1.32,Y,"PHOG0103256","P27875_BACSU"
6371893,P63275_HELPY,SwissProt,"Succinyl coA transferase",PHOG0273856...
6371893,P63275_HELPY,SwissProt,"Putative transferase",PHOG0273856,"Q4...
6371894,P36725_HELPY,GO_cellular_localization,"outer membrane",IDA,"P...
6371894,P36725_HELPY,GO_cellular_localization,"cytoplasm",N,"PHOG0372...
In this example, the numbers like 6371892 are JGI identifiers, which
were given by the user in the input FASTA file.

Note that in rows 4 & 5 of this example, there are two different
SwissProt descriptions of the same sequence 6371892 coming from the same
PHOG, PHOG0273856.  Instead of trying to compute a consensus description,
we can treat the description the same as the other annotations and simply
list each description in a separate row with all the PHOGs and orthologs
for which it occurs.
"""

from optparse import OptionParser
import sys

from pfacts003.phylofacts.models import UniProt, TreeNode,\
    OrthologTypes, preset_thresholds, UniProtEC, KEGG_Map_EC, GO_Term,\
    UniProtGO
from pfacts003.phog.orthologs import getOrthologQuerySet as\
    get_ortholog_query_results



def add_entry(annot_dict, annot_type, annot, evidence, phog, orth):
    """Add entry to dictionary
    
    Because defaultdict from the collections module is not available
    until python 2.5, this function is a work-around for default values
    (set, set).
    """
    key = (annot_type,annot,evidence)
    if not annot_dict.has_key(key):
        annot_dict[key] = (set(), set())
    annot_dict[key][0].add(phog)
    annot_dict[key][1].add(orth)


def print_entries(annot_dict, jgi, uniprot, ortholog_type, threshold):
    """Print annot_dictionary in expected CSV format
   
    This function walks through dictionary 'annot_dict', and print results
    in the requested CSV format as specified in the 'PHOG High-
    Throughput' Google Doc.

    This will be called several times as we print results for each
    input as we go (and thus can track progress by tailing the
    output of this program). Therefore, column headers are printed
    elsewhere - and not within this function.
    """
    for key,val in annot_dict.items():

        annot_type, annot, evidence = key
        phogs, orthologs = val
        
        # get_accession is a print function only, not a db query
        phogs = [p.get_accession(
                    ortholog_type=ortholog_type,
                    threshold=threshold)
                        for p in phogs]

        print ','.join([jgi,
                        uniprot,
                        annot_type,
                        annot,
                        evidence,
                        '"%s"' % (','.join(phogs)),
                        '"%s"' % (','.join(orthologs)),
                       ])

def generate_report(filename, ortholog_type, threshold):

    """Generate Report for input file

    This function is the heart of this module. Currently, there are some
    procedures that are used that less efficient for loading larger
    blocks of data. That is, an individual function call for a single
    amount of data is less efficient than if we can get all of that data
    in a single call.

    However, some of the other function that we have attempted, although
    more efficient, have proven to miss some data that is presented in
    the current web page.

    Additionally, the evidence column for KEGG data is not currently
    being retrieved.
    """
    phog_count = 0
    f = open(filename, 'r')
    for row in f:
        row = row.strip().split(',')
        if row[0] == 'JGI Accession':
            continue
        if row[0].endswith('sequences found in phogs'):
            continue
        jgi, uni = row

        (uniprot_object, phogs, best_phogs, orthologs,
         phog_ortho_dict, message) =\
            get_ortholog_query_results(uni, ortholog_type, threshold=threshold)
    
        entry_dict = {}

        # We want to know if a protein is found in an orthology group,
        # even if there are no annotations for that protein
        if phogs is not None and len(phogs) > 0:
            phog_count = phog_count + 1

        # There's no point checking further if no orthologs are returned
        if orthologs == None:
            continue

        for ortholog in orthologs:
            uniprot = ortholog.sequence_header.uniprot
            if uniprot is None:
                continue
            # get EC
            for ec in UniProtEC.objects.filter(
                uniprot=uniprot,
                # If in brenda, this is placed in the evidence column
                # commenting the following line
                #is_in_brenda_f=True,
            ):
                add_entry(entry_dict,
                      'EC',
                      '"%s: %s"' % (ec.ec, ec.ec.description),
                      ec.is_in_brenda_f and 'In Brenda' or 'Not in Brenda',
                      phog_ortho_dict[ortholog],
                      ec.uniprot.uniprot_identifier,
                )


            # Get KEGG

            # The 'get_kegg_map_ids' is an expensive way to get this as
            # we have a database query for each uniprot id. We also are
            # re-asking # for EC data that we already have. Leaving this
            # for now, but should # be updated. We also can't get
            # is_in_breda_data without expensive queries as this is
            # structured
            for kegg_map in uniprot.get_kegg_map_ids():
                add_entry(entry_dict,
                    'KEGG',
                    '"%s: %s"' % (kegg_map.id, kegg_map.title),
                    'TODO: is_in_brenda_f values',
                    phog_ortho_dict[ortholog],
                    uniprot.uniprot_identifier,
                )
        
            # get GO
            for u_g in UniProtGO.objects.filter(
                    uniprot=uniprot,
                    go_evidence__priority__lt=6,
                ):
                    for go in GO_Term.objects.filter(
                        term_type__in=(
                            'molecular_function',
                            'biological_process',
                            'cellular_component'
                        ),
                        uniprotgo=u_g,
                    ):
                        add_entry(entry_dict,
                            'GO_%s' % go.term_type,
                            '"%s %s"' % (go.acc, go.name),
                            u_g.go_evidence.evidence,
                            phog_ortho_dict[ortholog],
                            uniprot.uniprot_identifier,
                        )
        
       
            if uniprot.in_swissprot_f:
                add_entry(entry_dict,
                    'Swissprot',
                    '"%s"' % (uniprot.de.rstrip('; ')),
                    'in_swissprot_f',
                    phog_ortho_dict[ortholog],
                    uniprot.uniprot_identifier,
                )

        print_entries(entry_dict, jgi, uni, ortholog_type, threshold)

        # Delete dictionary as we may be processing very large files and
        # don't need to keep the results as they are now printed
        # (although possibly not yet flushed, to standard out)
        del entry_dict

    return phog_count


if __name__ == '__main__':

    usage = """usage: %prog [options] <input_filename>
    use --help option to see additional options
    """

    report_type = {'super': OrthologTypes.SuperOrtholog,
                   'tight': OrthologTypes.PHOG_T_Tight,
                   'medium': OrthologTypes.PHOG_T_Medium,
                   'loose': OrthologTypes.PHOG_T_Loose,
                   'custom': OrthologTypes.PHOG_T_Custom
                  }

    parser = OptionParser()
    parser.add_option(
        "-t", "--ortholog_type", dest="ortholog_type",
        default="tight",
        help="Choices: %s" % report_type.keys())
    parser.add_option(
        "-d", "--distance", dest="distance",
        default=None,
        help="Custom evolutionary distance (takes longer to run)")
     
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("%s\n" % usage)
    else:
        filename = args[0]

    if options.ortholog_type not in report_type.keys():
        parser.error("ortholog_type be one of the following: %s\n" %
                      report_type.keys())

    if options.distance != None and options.ortholog_type != 'custom':
        parser.error(("Custom distances can only be used with "
                      "--ortholog_type='custom' option"))

    if options.distance is not None:
        try:
            options.distance = float(options.distance)
        except ValueError:
            parser.error("Distance should be a floating point value. "
                         "Received '%s'" %options.distance)

    if options.ortholog_type == 'custom':
        options.threshold = options.distance
    else:
        options.threshold =\
                   preset_thresholds[report_type[options.ortholog_type]]

    sys.stderr.write(
        "Generating %s report with threshold %f\n" % (
            options.ortholog_type,
            options.threshold))


    print ("Input Identifier,UniProt Identifier,Annotation Type,"
           "Annotation,Experimental Evidence?,PHOGs,Orthologs")
    print ",,Total_in_PHOGs,,%d,," % generate_report(filename,
                  report_type[options.ortholog_type], options.threshold)
