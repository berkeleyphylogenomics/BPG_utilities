#!/bin/sh

if [ ! -e $1 ]; then
    echo "\nFile $1 does not exist or cannot be found!\n"
    exit
fi

blastall -p 'blastp' -i $1 -d /clusterfs/ohana/external/UniProt/current/protein -o 'family.blast.hits'
grep "^..|" family.blast.hits  | cut -f2 -d"|" > family.blast.accessions
fastacmd -i family.blast.accessions -d /clusterfs/ohana/external/UniProt/current/protein -o 'family.blast.hits.fasta'
