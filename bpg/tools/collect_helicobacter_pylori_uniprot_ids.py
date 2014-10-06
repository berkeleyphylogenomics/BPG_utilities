#!/usr/bin/env python

import os, sys
from pfacts003.utils.id_patterns import uniprot_identifier_re3

def main():
  f = open(
    '/clusterfs/ohana/sandbox/more_uniprot_ids_of_helicobacter_pylori_jgi_ids.csv')
  lines = f.readlines()
  f.close()
  uniprot_ids_of_jgi_id = {}
  for line in lines:
    if line[0:4] == 'lcl|':
      fields = line.strip().split(',')
      jgi_id = fields[0]
      if jgi_id not in uniprot_ids_of_jgi_id:
        uniprot_ids_of_jgi_id[jgi_id] = set()
      for field in fields[1:]:
        uniprot_ids_of_jgi_id[jgi_id].add(field.replace('"',''))
  for jgi_id in uniprot_ids_of_jgi_id:
    found_canonical_id = False
    ids_from_swissprot = set()
    ids_with_HELPY = set()
    for id in uniprot_ids_of_jgi_id[jgi_id]:
      accession, taxon = id.split('_')
      if taxon == 'HELPY':
        ids_with_HELPY.add(id)
      if uniprot_identifier_re3.match(id):
        ids_from_swissprot.add(id)
    ids_in_swissprot_and_HELPY = ids_from_swissprot & ids_with_HELPY
    if len(ids_in_swissprot_and_HELPY) > 0:
      print "%s,%s" % (jgi_id, list(ids_in_swissprot_and_HELPY)[0])
    elif len(ids_with_HELPY) > 0:
      print "%s,%s" % (jgi_id, list(ids_with_HELPY)[0])
    elif len(ids_from_swissprot) > 0:
      print "%s,%s" % (jgi_id, list(ids_from_swissprot)[0])
    else:
      print "%s,%s" % (jgi_id, list(uniprot_ids_of_jgi_id[jgi_id])[0])

if __name__ == '__main__':
  main()
