#!/bin/env python

"""
The download and insertion of Genomic Locus information was a project that was
started, and then stopped as other priorities pushed this off the plate.

Thus, the work here is incomplete. However, to preserve the work thus done,
this is checked into the repository 'as-is.' Care should be taken not to
execute this in production without more testing, debugging and documentation.
"""

import datetime
import cPickle # remove
import sys

import psycopg2

from pfacts003.phylofacts.models import GeneInfo, UniProtGeneID, UniProtTaxonomy
from pfacts003.utils.credentials import get_credentials


DB_NAME="pfacts003_sandbox"
DB_USER="bpg_user"

status_choices = ('NA', 'INFERRED', 'MODEL', 'PREDICTED', 'PROVISIONAL', 'REVIEWED', 'VALIDATED')

class Gene(object):
    def __init__(self):
        self.is_empty = True
        self.taxon = None
        self.gene_id = None
        self.status = None
        self.gn_accession = None
        self.start_position = None
        self.end_position = None
        self.orientation = None

    def clear(self):
        self.__init__()

    def _handle_conflict(self, gene_info, taxon, gene_id, status,
                         gn_accession, start_position, end_position,
                         orientation):
        """Handle input file conflicts

        It is possible for the input file to have duplicate start and
        end position records.

        One such example is:
            taxid   GeneID     status    start   end   orientation
            9606   100500719      -      89453  90011      -
            9606   100500719   INFERRED    100    658      +
        
        The above is a pseduo-gene and one we will not import.

        Another example is:
            taxid     GeneID  status     start      end    orientation
             3702    5007813 REVIEWED  21445902  21447340      -
             3702    5007813 REVIEWED  21445902  21447340      -
        
        Although the above appears to be a duplicate record, there are
        fields that we are not importing. Showing the other field
        entries in the same row order shows how the entries are
        different:
       
        taxid   GeneID   RNA_nucl_accession RNA_nucl_gi  protein_access.
        3702    5007813  NM_001084269.2      186491663    NP_001077738.1 
        3702    5007813  NM_001124039.1      186491666    NP_001117511.1
      
        """
        # Determine if records are identical:
        if self.gene_id == gene_id and\
           self.status == status and\
           self.start_position == start_position and\
           self.end_position == end_position and\
           self.orientation == orientation:
            # For our concerns, these records are identical
            # so, it doesn't matter
            return

        # If we don't have the gene record, we don't care about this
        # record anyhow
        conn = psycopg2.connect('dbname=%s host=db user=%s password=%s' % (
            DB_NAME, DB_USER, get_credentials(DB_USER)
        ))

        cur = conn.cursor()
        cur.execute("select uniprot_id from uniprot_gene_id where geneid=%d" % gene_id);
        uniprot = cur.fetchone()

        # If we do have the gene record, we need to take further steps
        if uniprot is not None:
            if self.gn_accession.startswith('NC_') and gn_accession.startswith('NW_'):
                # NC trumps NW
                return

            if self.gn_accession.startswith('NW_') and gn_accession.startswith('NC_'):
                # update with new and return
                self._overwrite(gene_info, taxon, gene_id, status,
                                gn_accession, start_position, end_position,
                                orientation)
                return

            print "Warning: Conflict in tax_on %s, gene %s. Start Position: %d fighting start position: %d" % (self.taxon, self.gene_id, start_position, self.start_position)
            #import pprint
            #pprint.pprint(gene_info[(taxon, gene_id)])


    def _overwrite(self, gene_info, taxon, gene_id, status, gn_accession, start_position, end_position, orientation):
            self.is_empty = False
            self.taxon = taxon
            self.gene_id = gene_id
            self.status = status
            self.gn_accession = gn_accession
            self.start_position = start_position
            self.end_position = end_position
            self.orientation = orientation

    def collate(self, gene_info, taxon, gene_id, status, gn_accession, start_position, end_position, orientation):
        # Quick Sanity checks
        if self.taxon is not None and self.taxon != taxon:
            raise Exception("Can't collate taxon (%s) with %s." % (taxon, self.taxon))

        if self.gene_id is not None and self.gene_id != gene_id:
            raise Exception("Can't collate gene (%s) with %s." % (gene_id, self.gene_id))

        if end_position is None and start_position is not None:
            raise Exception("End position needed if start position given")

        if self.start_position is not None and start_position is not None:
            # At this point, there is a conflict of locus information.
            # See see._handle_conflicts doc string
            self._handle_conflict(gene_info, taxon, gene_id, status, gn_accession, start_position, end_position, orientation)

        # Overwrite if more relevant data
        if self.start_position is None and start_position is not None:
            self._overwrite(gene_info, taxon, gene_id, status,
                            gn_accession, start_position, end_position,
                            orientation)

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (self.taxon, self.gene_id, self.status, self.start_position, self.end_position, self.orientation)

#def process_gene(gene, start_position,  end_position, orientation):
#    print u.uniprot, u.geneid, gene, start_position, end_position, orientation


def clean_field(field, type=str):
    field = field.strip()

    if field == '-':
        field = None

    if type == 'int' and field is not None:
        field = int(field)

    return field


def non_sorted_error(previous_number, current_number, tax_or_gene):

    if tax_or_gene == 'taxon':
        error_string = "Taxon id %d comes after %d." % (current_number,
                                                       previous_number)
    else:
        error_string = "Gene id %d comes after %d." % (current_number,
                                                       previous_number)

    sys.stderr.write("""
Error: The input file is not sorted by the first two columns as expected.
%s

At the time that this program was written, the input file was generated in
sorted order. To save memory and time, this program uses that fact and
streams through the data in sorted order. If the data is no longer sorted,
then this program  cannot give proper results.

Use the Unix shell sort program to sort numerically by the second column
and try again.

""" % error_string)
    sys.exit(1)


def parse_gene2acession(gene_info):
    f = open('/clusterfs/ohana/external/ncbi/gene/to_import/gene2accession', 'r')

    end_of_file = False
    previous_taxon = 0
    previous_gene = 0
    gene = Gene()

    # That is, while not end-of-file:
    while not end_of_file:
        line = f.readline()

        if line == "":
            end_of_file = True
            continue

        if line.startswith('#'):
            continue

        line = line.split('\t')

        taxon = clean_field(line[0], 'int')
        gene_id = clean_field(line[1], 'int')
        status = clean_field(line[2])
        gn_accession = clean_field(line[8])
        start_position = clean_field(line[9], 'int')
        end_position = clean_field(line[10], 'int')
        orientation = clean_field(line[11])

        # Sanity check that input file is sorted
        if taxon < previous_taxon:
            non_sorted_error(previous_taxon, taxon, tax_or_gene='taxon')
        elif taxon == previous_taxon and gene_id < previous_gene:
            non_sorted_error(previous_gene, gene_id, tax_or_gene='gene')

        # Collate (taxon, gene) information
        if taxon == previous_taxon and gene_id == previous_gene:
            gene.collate(gene_info, taxon, gene_id, status, gn_accession, start_position, end_position, orientation)
        else:
            #if not gene.is_empty:
            #    print gene
            gene.clear() 
            previous_gene = gene_id
            previous_taxon = taxon

    #    #if status is None and start_position is not None:
    #    #    print "Example demonstrates",gene, status, start_position, end_position
    #
    #    #print line
    #    #t = UniProtTaxonomy.objects.get(id=int(line[0]))
    #    #try:
    #    #    u = UniProtGeneID.objects.get(geneid=gene)
    #    #except:
    #    #    u = None
    #    ##print "Taxon: %d:%s" % (int(line[0]), t),
    #    #print "taxid: ",line[0],
    #    ##print "Gene ID: %d/UniProt id: %s" % (int(line[1]), u.uniprot),
    #    #print "gene", line[1],
    #    #print "status: ", line[2],
    #    #print "line[3]:", line[3]
    #    ##print line[0], t, gene, line[2], start_position, end_position, orientation
    #    #print
    #
    #    #process_gene(gene, start_position,  end_position, orientation)


def parse_gene_info():
    """Build fast-lookup for chromosome information

    The data for gene locus information is stored in the gene2accession
    file as downloaded from ncbi.  The data for which chromosome that
    gene is located is in the gene_info file. 

    We obviously need to join this data -- but only the data that we
    need. We only need a subset of this data to be placed within the
    database.

    Thus, it is fastest to sort through one file while reviewing the
    other. The fastest approach is an in-memory lookup. Of the two
    files, the gene_info is the easiest as there is less logic for
    resolving duplicate entries.

    So, first, the contents of the gene_info file is parsed into a
    lookup-dictionary based upon tuple keys: (taxon, gene_id).

    Although this input data is growing exponentially over time,
    currently this approach is fast and successful. This look-up table
    fits within 2-3 Gigs of resident memory at the time this module
    was written.  If the input file grows faster than memory is
    allocated to our node machines, a different approach will need to
    be taken.

    Currently, it takes approximately 7 minutes to create this
    dictionary, from the input file, on a node on the Ohana cluster.
    """
    f = open('/clusterfs/ohana/external/ncbi/gene/to_import/gene_info', 'r')

    end_of_file = False
    gene_info = {}

    print "Begin parsing gene_info: ", datetime.datetime.now()
    # That is, while not end-of-file:
    while not end_of_file:
        line = f.readline()

        if line == "":
            end_of_file = True
            continue

        if line.startswith('#'):
            continue

        line = line.split('\t')

        taxon = clean_field(line[0], 'int')
        gene_id = clean_field(line[1], 'int')
        locus_tag = clean_field(line[3])
        chromosome = clean_field(line[6])
        modification_date = clean_field(line[14])

        if (taxon, gene_id) in gene_info:
            print """
            Error: Input file no longer contains unique fields

            The method of using a tuple of (taxon, gene_id) as a unique key
            works because the input file did not have any duplicate entries for
            these fields.  This is no longer true. There was already a taxon
            and gene when this entry was attempted:

            Taxon: %s
            Gene: %s
            Locus Tag: %s
            Chromosome: %s
            Modification Date: %s

            This assumption is obviously no longer true and this part of the
            program needs rewritten.
            """ %(taxon, gene_id, locus_tag, chromosome, modification_date)
        gene_info[(taxon, gene_id)] = (locus_tag, chromosome, modification_date)

    f.close()
    print "Finished parsing gene_info: ", datetime.datetime.now()
    return gene_info

#parse_gene2acession(parse_gene_info())
parse_gene2acession({})
