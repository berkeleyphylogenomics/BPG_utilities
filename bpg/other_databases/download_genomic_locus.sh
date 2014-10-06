#!/bin/sh -l

# The download and insertion of Genomic Locus information was a project that was started, and
# then stopped as other priorities pushed this off the plate.
#
# Thus, the work here is incomplete. However, to preserve the work thus done, this is checked
# into the repository 'as-is.' Care should be taken not to execute this in production without
# more testing, debugging and documentation.
#

BLASTDB="/clusterfs/ohana/external"

cd $BLASTDB/tmp_downloads

wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz
wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
