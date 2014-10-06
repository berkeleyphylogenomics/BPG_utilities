#!/usr/bin/env python

import MySQLdb
from connect_db import *
from parse_orthomcl_genomes import get_taxon_of_abbrev

db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, \
                      user=phylofacts_db_user, passwd=phylofacts_db_passwd )
cur = db.cursor( MySQLdb.cursors.DictCursor )

def main():
  taxon_of_abbrev = get_taxon_of_abbrev()
  f = open("/usr/local/blastdb/orthomcl/groups_orthomcl-2.txt")
  lines = [line.rstrip() for line in f.readlines()]
  f.close()
  line_num = 0
  for line in lines:
    line_num += 1
    fields = line.split()
    orthology_group = fields[0].rstrip(':')
    for ortholog_spec in fields[1:]:
      abbrev, ortholog = ortholog_spec.split('|')
      taxon = taxon_of_abbrev[abbrev]
      sql = """REPLACE INTO orthomcl
                        SET orthology_group = '%s',
                            taxon_id = %s,
                            sequence_identifier = '%s'""" \
            % (orthology_group, taxon, ortholog)
      cur.execute(sql)
    if line_num % 1000 == 0:
      print "Processed %d lines" % line_num

if __name__ == '__main__':
  main()
