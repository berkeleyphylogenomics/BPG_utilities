#!/usr/bin/env python

import MySQLdb
from connect_db import *

def main():
  db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, \
                        user=phylofacts_db_user, passwd=phylofacts_db_passwd )
  cur        = db.cursor( MySQLdb.cursors.DictCursor )
  sql = """SELECT DISTINCT protein_sequence_id, 
                  chars 
            FROM  tree_node, 
                  protein_sequence, 
                  sequence 
            WHERE superorthologous_node_id <> 0 
              AND protein_sequence_id <> 0 
              AND protein_sequence_id = protein_sequence.id 
              AND sequence_id = sequence.id"""
  n_rows = cur.execute(sql)
  rows = cur.fetchall()
  for row in rows:
    print ">bpgprotseq%s" % row['protein_sequence_id']
    print row['chars']

if __name__ == '__main__':
  main()
