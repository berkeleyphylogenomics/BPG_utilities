#!/usr/bin/env python

## find_thresholded_phogs.py
## Author: Ruchira S. Datta
## Copyright (c) 2009, Regents of the University of California
## All rights reserved.
##
## Redistiribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
## o Redistributions of source code must retain the above copyright notice,
## this list of conditions and the following disclaimer.
##
## o Redistributions in binary form must reproduce the above copyright notice,
## this list of conditions and the following disclaimer in the documentation
## and/or other materials provided with the distribution.
##
## o Neither the name of the University of California, Berkeley nor the names
## of its contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY
## EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR 
## ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
## OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
## DAMAGE.

import MySQLdb
from connect_db import *
from optparse import OptionParser

def get_duplication_distance(cur, tree_node_id):
  sql = """SELECT MAX((IF(next_nearest_distance >= 0.0,
                          distance + next_nearest_distance,
                          10000)
                      )/2)
              AS  duplication_distance
            FROM  maximal_tree_node_nearest_protein_sequence_with_species
            WHERE tree_node_id = %s
        """ % tree_node_id
  n_rows = cur.execute(sql)
  if n_rows == 0:
    return -1.0
  return float(cur.fetchone()['duplication_distance'])

def get_difference_distance(cur, tree_node_id):
  sql = """SELECT MAX((IF(next_nearest_distance >= 0.0,
                          next_nearest_distance - distance,
                          10000)
                      )/2)
              AS  difference_distance
            FROM  maximal_tree_node_nearest_protein_sequence_with_species
            WHERE tree_node_id = %s
        """ % tree_node_id
  n_rows = cur.execute(sql)
  if n_rows == 0:
    return -1.0
  return float(cur.fetchone()['difference_distance'])

class ProximalTreeNode:
  def __init__(self, cur, visit_id, tree_node_id, left_id, right_id, 
                parent_visit_id, verbose=False):
    if visit_id > 0:
      duplication_distance = get_duplication_distance(cur, tree_node_id)
    else:
      duplication_distance = 5000.0
    if visit_id > 0:
      difference_distance = get_difference_distance(cur, tree_node_id)
    else:
      difference_distance = 10000.0
    if verbose:
      print "New ProximalTreeNode: "
      print "  visit_id: %d, tree_node_id: %s, left_id: %s, right_id: %s, " \
        % (visit_id, tree_node_id, left_id, right_id),
      print "parent_visit_id: %d, duplication_distance: %g, " \
        % (parent_visit_id, duplication_distance),
      print "difference distance: %g, " % difference_distance
    self.visit_id = visit_id
    self.tree_node_id = tree_node_id
    self.left_id = left_id
    self.right_id = right_id
    self.parent_visit_id = parent_visit_id
    self.children_visit_ids = set()
    self.duplication_distance = duplication_distance
    self.difference_distance = difference_distance
    self.has_maximal_descendant_at_threshold = {}
    self.has_maximal_ancestor_at_threshold = {}
    self.is_phogt_at_threshold = {}

def get_tree_id(cur, bpg_accession, method='nj'):
  sql = """SELECT tree.id AS tree_id
            FROM  book,
                  alignment,
                  tree
            WHERE book.id = book_id
            AND   alignment.id = alignment_id
            AND   scopid_or_bpgid = '%s'
            AND   method = '%s'
        """ % (bpg_accession, method)
  n_rows = cur.execute(sql)
  if n_rows == 0:
    return None
  else:
    return cur.fetchone()['tree_id']

def print_all_proximal_nodes(cur, tree_id):
  sql = """SELECT     tree_node_id,
                      left_id,
                      right_id
            FROM      tree_node,
                      maximal_tree_node_nearest_protein_sequence_with_species
            WHERE     tree_node.id = tree_node_id
            AND       tree_id = %s
            GROUP BY  left_id
        """ % tree_id
  n_rows = cur.execute(sql)
  print "%d proximal nodes" % n_rows
  print "-------------------------------------"
  print "| tree_node_id | left_id | right_id |"
  print "-------------------------------------"
  for row in cur.fetchall():
    print "|    %09d | %07d |  %07d |" \
      % (row['tree_node_id'], row['left_id'], row['right_id'])
    print "-------------------------------------"
  print
    
def get_breadth_first_proximal_nodes(cur, tree_id, verbose=False):
  breadth_first_search_order = {}
  current_index = 0
  num_nodes = 0
  fake_root = ProximalTreeNode(cur, num_nodes, -1, -1, -1, -1, verbose)
  breadth_first_search_order[0] = fake_root
  num_nodes += 1
  sql = """SELECT     tree_node_id,
                      left_id,
                      right_id
            FROM      tree_node,
                      maximal_tree_node_nearest_protein_sequence_with_species
            WHERE     tree_node.id = tree_node_id
            AND       tree_id = %s
            ORDER BY  left_id
            LIMIT     1
        """ % tree_id
  n_rows = cur.execute(sql)
  while n_rows > 0:
    row = cur.fetchone()
    breadth_first_search_order[num_nodes] \
      = ProximalTreeNode(cur, num_nodes, row['tree_node_id'],
                          row['left_id'], row['right_id'], 0, verbose)
    fake_root.children_visit_ids.add(num_nodes)
    num_nodes += 1
    sql = """SELECT     tree_node_id,
                        left_id,
                        right_id
              FROM      tree_node,
                        maximal_tree_node_nearest_protein_sequence_with_species
              WHERE     tree_node.id = tree_node_id
              AND       tree_id = %s
              AND       left_id > %s
              ORDER BY  left_id
              LIMIT     1
          """ % (tree_id, row['right_id'])
    n_rows = cur.execute(sql)
  current_index += 1
  while current_index < num_nodes:
    current_node = breadth_first_search_order[current_index]
    sql = """SELECT     tree_node_id,
                        left_id,
                        right_id
              FROM      tree_node,
                        maximal_tree_node_nearest_protein_sequence_with_species
              WHERE     tree_node.id = tree_node_id
              AND       tree_id = %s
              AND       left_id > %s
              AND       right_id < %s
              ORDER BY  left_id
              LIMIT     1
          """ % (tree_id, current_node.left_id, current_node.right_id)
    n_rows = cur.execute(sql)
    while n_rows > 0:
      row = cur.fetchone()
      breadth_first_search_order[num_nodes] \
        = ProximalTreeNode(cur, num_nodes, row['tree_node_id'],
                            row['left_id'], row['right_id'],
                            current_index, verbose)
      current_node.children_visit_ids.add(num_nodes)
      num_nodes += 1
      sql = """SELECT     tree_node_id,
                          left_id,
                          right_id
                FROM      tree_node,
                      maximal_tree_node_nearest_protein_sequence_with_species
                WHERE     tree_node.id = tree_node_id
                AND       tree_id = %s
                AND       left_id > %s
                AND       right_id < %s
                ORDER BY  left_id
                LIMIT     1
            """ % (tree_id, row['right_id'], current_node.right_id)
      n_rows = cur.execute(sql)
    current_index += 1
  if verbose:
    for i in xrange(num_nodes):
      print "Parent: %d " % i,
      print "Children: ", breadth_first_search_order[i].children_visit_ids
      print " Duplication distance: ",
      print breadth_first_search_order[i].duplication_distance
      print " Difference distance: ",
      print breadth_first_search_order[i].difference_distance
    print "Put %d proximal nodes (incl. 1 fake root) in order" % num_nodes
  return (breadth_first_search_order, num_nodes)
              
def make_distances_monotonic(cur, breadth_first_search_order, num_nodes,
                              verbose = False):
  for i in xrange(num_nodes):
    j = num_nodes - 1 - i
    parent_visit_id =  breadth_first_search_order[j].parent_visit_id
    if parent_visit_id > 0:
      if breadth_first_search_order[parent_visit_id].duplication_distance \
          < breadth_first_search_order[j].duplication_distance:
        if verbose:
          print "Changing duplication distance of %d " % parent_visit_id,
          print " from %g to %g " \
                % (breadth_first_search_order[
                                        parent_visit_id].duplication_distance,
                    breadth_first_search_order[j].duplication_distance)
        breadth_first_search_order[parent_visit_id].duplication_distance \
          = breadth_first_search_order[j].duplication_distance

def put_distances_in_tree_node_table(cur, breadth_first_search_order,
                                    num_nodes, verbose = False):
  # Skip the fake root node at index 0
  for i in xrange(1,num_nodes):
    node = breadth_first_search_order[i]
    sql = """UPDATE tree_node
                SET duplication_distance = %g
              WHERE id = %s
          """ % (node.duplication_distance, node.tree_node_id)
    if verbose:
      print sql
    n_rows = cur.execute(sql)
    sql = """UPDATE tree_node
                SET difference_distance = %g
              WHERE id = %s
          """ % (node.difference_distance, node.tree_node_id)
    if verbose:
      print sql
    n_rows = cur.execute(sql)

def compute_phogts_at_threshold(breadth_first_search_order, num_nodes, 
                                threshold, kind, verbose):
  if verbose:
    print "Finding PHOG-T's at %s threshold %g" % (kind, threshold)
  for i in xrange(num_nodes):
    node = breadth_first_search_order[i]
    node.has_maximal_descendant_at_threshold[kind] = False
    node.has_maximal_ancestor_at_threshold[kind] = False
    node.is_phogt_at_threshold[kind] = False
  phogts_of_kind = set()
  # skip the fake root at 0
  for i in xrange(num_nodes-1):
    j = num_nodes - 1 - i
    node = breadth_first_search_order[j]
    parent = breadth_first_search_order[node.parent_visit_id]
    if node.has_maximal_descendant_at_threshold[kind]:
      parent.has_maximal_descendant_at_threshold[kind] = True
    else:
      if node.duplication_distance >= threshold:
        if verbose:
          print "%s PHOG-T at %s (%s, %s) %g" % (kind, node.tree_node_id,
                                              node.left_id, node.right_id,
                                              node.duplication_distance)
        phogts_of_kind.add(j)
        node.is_phogt_at_threshold[kind] = True
        parent.has_maximal_descendant_at_threshold[kind] = True
  return phogts_of_kind

def put_phogts_in_tree_node_table(cur, tree_id, breadth_first_search_order, 
                            num_nodes, phogts, threshold, kind, verbose = False):
  if breadth_first_search_order[1].left_id <> 1:
    sql = """SELECT id
               FROM tree_node
              WHERE tree_id = %s
              AND   left_id = 1
          """ % tree_id
    n_rows = cur.execute(sql)
    root_node_id = cur.fetchone()['id']
    if verbose:
      print "Root node: ", root_node_id
    sql = """UPDATE tree_node
                SET phogt_%s_node_id = %s
              WHERE tree_id = %s
              AND   protein_sequence_id <> 0
          """ % (kind, root_node_id, tree_id)
    if verbose:
      print sql
    cur.execute(sql)
  for i in xrange(1, num_nodes):
    node = breadth_first_search_order[i]
    if node.duplication_distance >= threshold:
      sql = """UPDATE tree_node
                  SET phogt_%s_node_id = %s
                WHERE tree_id = %s
                AND   left_id >= %d
                AND   right_id <= %d
                AND   protein_sequence_id <> 0
            """ % (kind, node.tree_node_id, tree_id, 
                  node.left_id, node.right_id)
      if verbose:
        print sql
      n_rows = cur.execute(sql)
  sql = """UPDATE tree_node
              SET is_phogt_%s = FALSE
            WHERE tree_id = %s
        """ % (kind, tree_id)
  n_rows = cur.execute(sql)
  tree_node_ids = [breadth_first_search_order[i].tree_node_id for i in phogts]
  sql = """UPDATE tree_node
              SET is_phogt_%s = TRUE
            WHERE id IN (%s)
        """ % (kind, ','.join(['%d' % id for id in tree_node_ids]))
  if verbose:
    print sql
  n_rows = cur.execute(sql)

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-m", "--method", dest="method", default="nj",
                      help="tree method for which to find PHOG-T's.")
  opt_parser.add_option("-v", "--verbose", dest="verbose",
              action="store_true", default=True,
              help="Whether to print verbose output")
  opt_parser.add_option("-q", "--quiet", dest="verbose",
              action="store_false", default=True,
              help="Whether to suppress verbose output")
  (options, args) = opt_parser.parse_args()
  if len(args) < 1:
    opt_parser.error("Incorrect number of arguments")
  else:
    bpg_accession = args[0]
  verbose = options.verbose
  db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, \
                        user=phylofacts_db_user, passwd=phylofacts_db_passwd, \
                        charset='utf8' )
  cur = db.cursor( MySQLdb.cursors.DictCursor )
  tree_id = get_tree_id(cur, bpg_accession, options.method)
  if not tree_id:
    opt_parser.error('Unable to find %s tree for book %s in database'
                      % (method, bpg_accession))
  if verbose:
    print_all_proximal_nodes(cur, tree_id)
  (breadth_first_search_order, num_proximal_nodes) \
    = get_breadth_first_proximal_nodes(cur, tree_id, verbose)
  make_distances_monotonic(cur, breadth_first_search_order, num_proximal_nodes,
                            verbose)
  put_distances_in_tree_node_table(cur, breadth_first_search_order, 
                            num_proximal_nodes, verbose)
  thresholds = {}
  thresholds['tight'] = 0.09375
  thresholds['medium'] = 0.296874
  thresholds['loose'] = 0.9375
  for kind in thresholds.keys():
    phogts = compute_phogts_at_threshold(breadth_first_search_order,
                          num_proximal_nodes, thresholds[kind], kind, verbose)
    put_phogts_in_tree_node_table(cur, tree_id, breadth_first_search_order, 
                                  num_proximal_nodes, phogts, thresholds[kind],
                                  kind, verbose)

if __name__ == '__main__':
  main()
