#!/bin/sh -l

# Download and format UniProt - uniprot_sprot and uniprot_trembl

# Do work in working directory.

BLASTDB="/clusterfs/ohana/external"
cd $BLASTDB/tmp_downloads

log=$(pwd)/go_download.log

echo "----------------------" > $log
date >> $log

# uniprot_sprot.  Format here is ftp://user:password@address

month=`date +"%Y%m"`

url="http://archive.geneontology.org/latest-full/go_$month-assocdb-tables.tar.gz"
file="assocdb-tables.tar.gz"
curl -s $url -o $file >> $log 2>&1

# Move to subdirectory - name is date.

today=`date +"%Y-%m-%d"`
dir=$BLASTDB/go/$today
mkdir -p $dir

mv $file $dir

# Relink current to new directory.  Need to delete old link first.

unlink $BLASTDB/go/current
ln -s $dir $BLASTDB/go/current

# Delete versions more than 3 months old.

find  $BLASTDB/go -type d -mtime +90 -print0 | xargs -0 rm -r -f >> $log 2>&1

# Untar the table files

cd $dir
tar xzvf $file >> $log 2>&1

# Remove the tar archive
rm $file

# Send mail and exit.

echo "update current" | mail -s "check Gene Ontology download on ohana" bpg.shailen@gmail.com

