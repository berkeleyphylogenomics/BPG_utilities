#!/bin/sh 


log=/clusterfs/ohana/external/tmp_downloads/pdb_structures_download.log


/usr/bin/rsync -rlpt -v -z --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ /clusterfs/ohana/external/pdb_structures/  > $log 2>&1 

mutt bpg.shailen@gmail.com -s "PDB Structures Rsync'ed" -a $log << EOT
The PDB Structures were just rsync'ed. 

Attached is a copy of the rsync output/log file:
* Normal output log ($log)

Please review the last entries of these logs for any errors. If structures did
change this week, these need to be manually updated on Limahuli (until that
system is phased out). The Limahuli system is not set up in such an easy way to
update.


EOT
