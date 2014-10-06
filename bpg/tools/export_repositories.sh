#!/bin/bash
# Small Bash Script to take code export archive

REPOSITORIES="bpg_phylofacts bpg_research bpg_static bpg_teaching"

for repo in $REPOSITORIES; do
   svn export https://bpg.unfuddle.com/svn/$repo
   echo "Compressing $repo"
   tar -cvjf $repo.tar.bz2 $repo
   echo "Removing $repo"
   rm -Rf $repo
done
