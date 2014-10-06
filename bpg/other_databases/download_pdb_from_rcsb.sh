#!/bin/sh -l

# Download and format pdb db sequences from rcsb.

# Do work in working directory.
cd $BLASTDB/tmp_downloads

log=pdb_rcsb_download.log

echo "----------------------" >> $log
date >> $log

# Download single zipped file into this directory
url="ftp.wwpdb.org/pub/pdb/derived_data"
curl -s $url/pdb_seqres.txt.gz -o pdb_rcsb_full.gz >> $log  2>&1

# Unzip into this directory.
echo Unzipping... >> $log 2>&1
gunzip -f pdb_rcsb_full.gz >> $log 2>&1

# Create two additional, smaller DBs: one without DNA (pdb_rcsb),
# and one with identical seqs filtered out (pdb_rcsb_nr)
filter_out_dna.py pdb_rcsb_full > pdb_rcsb
uniqueseq pdb_rcsb_nr -db pdb_rcsb >> $log 2>&1
# Modify name given by uniqueseq.
mv pdb_rcsb_nr.seq pdb_rcsb_nr

# Format each of the three DBs.
#   -i input
#   -o T (for "True"): parse SeqId and create indices.
#   -l logfile

formatdb -i pdb_rcsb_full -o T -l pdb_rcsb_full_formatdb.log -p T >> $log 2>&1
formatdb -i pdb_rcsb -o T -l pdb_rcsb_formatdb.log -p T >> $log 2>&1
formatdb -i pdb_rcsb_nr -o T -l pdb_rcsb_nr_formatdb.log -p T >> $log 2>&1

# If seems to have worked (pdb_rcsb_nr.pin exists)
# move previous files out of way
# and move these files to BLASTDB directory.

file_list="pdb_rcsb pdb_rcsb_full.phr pdb_rcsb_nr pdb_rcsb_nr.psi pdb_rcsb.psi pdb_rcsb_full.pin pdb_rcsb_nr.psq  pdb_rcsb.psq pdb_rcsb_full.psd  pdb_rcsb_nr.phr pdb_rcsb.phr pdb_rcsb_full pdb_rcsb_full.psi  pdb_rcsb_nr.pin pdb_rcsb.pin pdb_rcsb_full.psq  pdb_rcsb_nr.psd pdb_rcsb.psd"
#pdb_rcsb_download.log  pdb_rcsb_nr_formatdb.log
#pdb_rcsb_formatdb.log pdb_rcsb_full_formatdb.log

if [ -e pdb_rcsb_nr.pin ]
then
# Move last download to <old> directory.
#   mv $BLASTDB/pdbaa* $BLASTDB/old
#   mv $BASE/pdb_rcsb_50nr* $BASE/test/old
  
    mv $BLASTDB/pdb_sequences/pdb_rcsb* $BLASTDB/old

# Move current download to <active> directory.
#   mv pdbaa pdbaa.* $BLASTDB
#   mv pdb_rcsb_50nr.* $BASE/test
    mv pdb_rcsb* $BLASTDB/pdb_sequences
fi

export REPLYTO='bpg.shailen@gmail.com'

mutt bpg.shailen@gmail.com@gmail.com -s "PDB Sequences Download/Installation" -a $log -a $log << EOT
The download_pdb_from_rcsb.sh script has just finished executing.

Attached is a copy of the log file:
* Normal output log ($log)

Please review the last entries of these logs for any errors.

EOT
