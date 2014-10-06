#!/usr/bin/env python

import os
import psycopg2
import cPickle
import psycopg2.extras
from pfacts003.utils.credentials import get_credentials

def main():
  dir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5'
  os.chdir(dir)
  info_of_uniprot_accession = {}
  f = open("all_uniprot_accessions.txt")
  for line in f.readlines():
    taxon_id, accession = line.strip().split(',')
    info_of_uniprot_accession[accession] = {}
    info_of_uniprot_accession[accession]['taxon'] = taxon_id
  f.close()

  uniprot_accessions_of_uniprot_id = {}

  bpg_password = get_credentials('bpg_user')

  connection = psycopg2.connect(
    "dbname='%s' user='%s' host='db' password='%s'" %
    ('pfacts003_test', 'bpg_user', bpg_password))

  cur = connection.cursor('server_side_cursor', 
                          cursor_factory = psycopg2.extras.DictCursor)

  sql = 'SELECT * FROM uniprot_dat_index'
  cur.execute(sql)
  for row in cur:
    if row[1] in info_of_uniprot_accession:
      accession = row[1]
      uniprot_id = int(row[2])
      info_of_uniprot_accession[accession]['uniprot_id'] \
          = uniprot_id
      if uniprot_id not in uniprot_accessions_of_uniprot_id:
        uniprot_accessions_of_uniprot_id[uniprot_id] = set()
      uniprot_accessions_of_uniprot_id[uniprot_id].add(accession)

  f = open("info_of_uniprot_accession.pkl", "w")
  cPickle.dump(info_of_uniprot_accession, f)
  f.close()
  f = open("uniprot_accessions_of_uniprot_id.pkl", "w")
  cPickle.dump(uniprot_accessions_of_uniprot_id, f)
  f.close()

if __name__ == '__main__':
  main()
