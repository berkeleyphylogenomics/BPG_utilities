#!/usr/bin/env python

import os

def main():
  dir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5'
  os.chdir(dir)
  f = open("sequences_without_uniprot_accessions.id")
  ensembl_ids = [line.strip().split(':')[1] for line in f.readlines()]
  uniprot_accession_of_ensembl_id = {}
  f.close()
  for ensembl_id in ensembl_ids:
    os.chdir('blast_results')
    os.chdir(ensembl_id)
    f = open('%s.blast_results' % ensembl_id)
    lines = f.readlines()
    f.close()
    for line in lines:
      if len(line) > 0 and line[0] == '#':
        continue
      break
    fields = line.strip().split()
    pwid = float(fields[2])
    if pwid >= 98.0:
      uniprot_accession_of_ensembl_id[fields[0].split(':')[1]] \
        = fields[1].split('|')[1]
    os.chdir('../..')

  print "EnsemblId,UniProtAccession"
  for ensembl_id in uniprot_accession_of_ensembl_id:
    print "%s,%s" % (ensembl_id, uniprot_accession_of_ensembl_id[ensembl_id])

if __name__ == '__main__':
  main()
