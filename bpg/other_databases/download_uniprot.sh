#!/bin/sh -l

# Download and format UniProt - uniprot_sprot and uniprot_trembl

# Do work in working directory.

BLASTDB="/clusterfs/ohana/external"
cd $BLASTDB/tmp_downloads

log=uniprot_download.log

echo "----------------------" > $log
date >> $log

# uniprot_sprot.  Format here is ftp://user:password@address

ftp="ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"
sprot="uniprot_sprot.fasta"
curl $ftp/$sprot.gz -o $sprot.gz >> $log 2>&1

sprot_dat="uniprot_sprot.dat"
curl $ftp/$sprot_dat.gz -o $sprot_dat.gz >> $log 2>&1

# uniprot_trembl
trembl="uniprot_trembl.fasta"
curl $ftp/$trembl.gz -o $trembl.gz >> $log 2>&1

trembl_dat="uniprot_trembl.dat"
curl $ftp/$trembl_dat.gz -o $trembl_dat.gz >> $log 2>&1

# keyword list
keywlist="keywlist.txt"
curl $ftp/docs/$keywlist -o $keywlist >> $log 2>&1

echo `date` >> $log
echo "Download complete" >> $log

# Unzip both, cat together into "sprot_trembl".

gunzip $sprot >> $log 2>&1
gunzip $trembl >> $log 2>&1
cat $sprot $trembl > protein
rm $sprot $trembl

gunzip $sprot_dat >> $log 2>&1
gunzip $trembl_dat >> $log 2>&1
cat $sprot_dat $trembl_dat > uniprot.dat
rm $sprot_dat $trembl_dat

# Format.
#   -i input
#   -o T (for "True"): parse SeqId and create indices.
#   -l logfile

formatdb -i protein -o T -l uniprot_format.log

# Move to subdirectory - name is date.

today=`date +"%Y-%m-%d"`
dir=$BLASTDB/UniProt/$today
mkdir -p $dir

mv protein* $dir
mv uniprot.dat $dir
mv keywlist.txt $dir

# Send mail and exit.

echo "update current" | mail -s "check UniProt download" bpg.shailen@gmail.com

# Relink current to new directory.  Need to delete old link first.

rm $BLASTDB/UniProt/to_import
ln -s $dir $BLASTDB/UniProt/to_import

# Delete versions more than 3 months old.

find  $BLASTDB/UniProt -type d -mtime +90 -print0 | xargs -0 rm -r -f >> $log 2>&1
