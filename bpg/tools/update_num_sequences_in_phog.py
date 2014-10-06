#!/usr/bin/env python

import MySQLdb
from connect_db import *
db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, \
                      user=phylofacts_db_user, passwd=phylofacts_db_passwd )
cur = db.cursor( MySQLdb.cursors.DictCursor )
select_cur = db.cursor( MySQLdb.cursors.DictCursor )
update_cur = db.cursor( MySQLdb.cursors.DictCursor )

sql = """SELECT id, tree_id, left_id, right_id, num_sequences
          FROM  tree_node 
          WHERE is_superorthologous = TRUE"""
print sql
n_phogs = cur.execute(sql)
print "Selected %d PHOGs" % n_phogs
num_updated_phogs = 0
for i in xrange(n_phogs):
  num_updated_phogs = num_updated_phogs + 1
  row = cur.fetchone()
  cur_num_seqs = row['num_sequences']
  sql = """SELECT COUNT(DISTINCT tree_node.id)
                  AS true_num_sequences
            FROM  tree_node
            WHERE tree_id = %s
            AND   left_id > %s
            AND   right_id < %s
            AND   protein_sequence_id <> 0
        """ % (row['tree_id'], row['left_id'], row['right_id'])
  print sql
  n_rows = select_cur.execute(sql)
  if n_rows > 0:
    num_seqs_row = select_cur.fetchone()
    true_num_seqs = num_seqs_row['true_num_sequences']
    if true_num_seqs != cur_num_seqs:
      sql = """UPDATE tree_node
                  SET num_sequences = %d
                WHERE id = %s
            """ % (true_num_seqs, row['id'])
      print sql
      update_cur.execute(sql)
