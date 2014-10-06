#!/bin/sh -l

# Download and format UniProt - uniprot_sprot and uniprot_trembl

# Do work in working directory.

BLASTDB="/clusterfs/ohana/external"
cd $BLASTDB/tmp_downloads

log=uniprot_taxonomy_download.log

echo "----------------------" > $log
date >> $log

# Download TSV file from UniProt

url="http://www.uniprot.org/taxonomy/?query=*&force=yes&format=tab"
table="taxonomy-all.tab"
curl -s $url -o $table >> $log 2>&1

# Move to subdirectory - name is date.

today=`date +"%Y-%m-%d"`
dir=$BLASTDB/taxonomy/$today
mkdir -p $dir

mv $table $dir

# Send mail and exit.

echo "update current" | mail -s "check UniProt taxonomy download on ohana" bpg.shailen@gmail.com

# Relink to_import to new directory.
ln -s $dir $BLASTDB/taxonomy/to_import

# Delete versions more than 3 months old.

find  $BLASTDB/taxonomy -type d -mtime +90 -print0 | xargs -0 rm -r -f >> $log 2>&1
