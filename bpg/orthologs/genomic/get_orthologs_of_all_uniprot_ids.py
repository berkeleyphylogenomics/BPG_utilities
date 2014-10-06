#!/usr/bin/env python

from mpi4py import MPI
import cPickle
import psycopg2, psycopg2.extras
from pfacts003.utils.credentials import get_credentials
import numpy
import os
import sys

( TAG_DATABASE_ROW,
  TAG_DATABASE_DONE,
  TAG_TREES_COMPUTED_ORTHOLOGS,
  TAG_ORTHOLOG_REQUEST,
  TAG_ORTHOLOG_RESPONSE,
  TAG_ORTHOLOG_DONE,
) = range(6)

class ServedTreeNode:
  def __init__(self, row):
    self.tree_id = int(row[0])
    self.left_id = int(row[1])
    self.right_id = int(row[2])
    self.lastIncomparableNode = None
    self.lastDescendantNode = None
    if numpy.isnan(row[3]):
      self.duplication_distance = None
    else:
      self.duplication_distance = row[3]
    if numpy.isnan(row[4]):
      self.greatest_duplication_distance_of_maximal_descendant = None
    else:
      self.greatest_duplication_distance_of_maximal_descendant = row[4]
    if numpy.isnan(row[5]):
      self.uniprot_id = None
    else:
      self.uniprot_id = int(row[5])
    self.parent = None
    self.children = set()

  def findRoot(self):
    if self.parent:
      return self.parent.findRoot()
    else:
      return self

  def countNodes(self):
    sum = 1
    for child in self.children:
      sum += child.countNodes()
    return sum

  # Method to attempt to insert a new node in the tree starting at self
  # This method returns whether or not the new node was inserted
  # This method should *not* be called after sortChildren, i.e., sortChildren
  # should only be called after all nodes have been inserted.
  def insertNode(self, node):
    if self.parent and (node.left_id < self.parent.left_id or \
          node.right_id > self.parent.right_id):
        return self.parent.insertNode(node)
    if node.left_id < self.left_id and node.right_id > self.right_id:
      # node is an ancestor of self
      if self.parent:
        # Now we know node.left_id >= self.parent.left_id and node.right_id <=
        # self.parent.left_id, so the node should be a descendant of the parent
        # Let p -> c mean p is the parent of c
        # Current: self.parent -> self
        # New: self.parent -> node -> self
        self.parent.children.remove(self)
        self.parent.children.add(node)
        node.parent = self.parent
      node.children.add(self)
      self.parent = node
      return True
    elif node.left_id > self.left_id and node.right_id < self.right_id:
      # node is a descendant of self
      self.lastDescendantNode = node
      child_list = list(self.children)
      isNewLoneSiblingOfCurrentChildren = True
      # self.children may change during the loop
      # specifically, node may be added to self.children, and multiple current
      # children of self may become children of node instead
      # on the other hand, node may be a descendant of one of self.children
      # instead, in which case self.children will not change
      for i in xrange(len(child_list)):
        child = child_list[i]
        # avoid an infinite loop: if the child found node was incomparable and
        # asked its parent to try to insert it instead, don't ask the child to
        # insert the node again
        if (child.lastIncomparableNode is None or \
            node.left_id != child.lastIncomparableNode.left_id) \
            and child.insertNode(node):
          isNewLoneSiblingOfCurrentChildren = False
      if isNewLoneSiblingOfCurrentChildren:
        # node was not inserted by any of the children
        node.parent = self
        self.children.add(node)
      return True
    else:
      # The node is neither an ancestor nor a descendant of self
      self.lastIncomparableNode = node
      if self.parent and (self.parent.lastDescendantNode is None or \
                          not self.parent.lastDescendantNode == node):
        return self.parent.insertNode(node)
      else:
        return False

  def pruneEmptyChildren(self):
    childrenToPrune = set()
    for child in self.children:
      child.pruneEmptyChildren()
    for child in self.children:
      if len(child.children) == 0 and child.uniprot_id is None:
        childrenToPrune.add(child)
    self.children = self.children - childrenToPrune

  # Method to sort children in order of left_id
  # The insertNode method should *not* be called after sortChildren, i.e.,
  # sortChildren should only be called after all nodes have been inserted.
  def sortChildren(self):
    child_list = list(self.children)
    child_list.sort(key = lambda child : child.left_id)
    self.children = child_list
    for child in self.children:
      child.sortChildren()

  def fillUniProtIds(self, uniprot_id_array, index_of_uniprot_id, node_of_index,
                    i, index_of_left_id_in_uniprot_id_array,
                    left_right_id_array, right_id_of_left_id, 
                    node_of_left_id, k):
    node_of_left_id[self.left_id] = self
    left_right_id_array[k] = self.left_id
    right_id_of_left_id[self.left_id] = self.right_id
    k += 1
    if self.uniprot_id is not None:
      if len(self.children) > 0:
        print "Weirdness! Node %d in tree %d has children and uniprot_id" \
          % (self.left_id, self.tree_id)
      uniprot_id_array[i] = self.uniprot_id
      index_of_uniprot_id[self.uniprot_id] = i
      node_of_index[i] = self
      index_of_left_id_in_uniprot_id_array[self.left_id] = i
      i += 1
    for child in self.children:
      (i, k) = child.fillUniProtIds(uniprot_id_array, index_of_uniprot_id, 
                              node_of_index, i,
                              index_of_left_id_in_uniprot_id_array,
                              left_right_id_array, right_id_of_left_id, 
                              node_of_left_id, k)
    left_right_id_array[k] = self.right_id
    k += 1
    return (i, k)

  def printSelf(self, level=0, outf=sys.stdout):
    for i in range(level):
      outf.write(' ')
    outf.write('Tree_id: %d Left_id: %d ' % (self.tree_id, self.left_id))
    if self.uniprot_id:
      outf.write('Uniprot_id: %d ' % self.uniprot_id)
    if self.greatest_duplication_distance_of_maximal_descendant:
      outf.write('Thresholds ')
      if self.greatest_duplication_distance_of_maximal_descendant == -1.0:
        outf.write('[0.0,')
      else:
        outf.write('(%f,' 
                  % self.greatest_duplication_distance_of_maximal_descendant)
      outf.write('%f] {\n' % self.duplication_distance)
      for child in self.children:
        child.printSelf(level=level+1, outf=outf)
      for i in range(level):
        outf.write(' ')
      outf.write('} Right_id: %d\n' % self.right_id)
    else:
      outf.write('\n')
      
class ServedTree:
  def __init__(self):
    self.root = None
    self.lastInsertedNodes = set()
    self.nodes_of_uniprot_id = {}
    self.total_num_inserted_nodes = 0
    self.uniprot_id_array = None
    self.thresholds_of_orthology = None
    self.phogs_supporting_orthology = None
    self.index_of_uniprot_id = {}
    self.num_uniprot_ids = 0
    self.end_uniprot_index = 0

  def addNodeFromRow(self, row):
    node = ServedTreeNode(row)
    self.total_num_inserted_nodes += 1
    if node.uniprot_id is not None:
      try:
        self.nodes_of_uniprot_id[node.uniprot_id].add(node)
      except KeyError:
        self.nodes_of_uniprot_id[node.uniprot_id] = set([node])
      self.num_uniprot_ids += 1
    if node.left_id == 1:
      self.root = node
      for lastInsertedNode in self.lastInsertedNodes:
        assert node.insertNode(lastInsertedNode)
      self.lastInsertedNodes = set([node])
    else:
      nodesWhereInserted = set()
      for lastInsertedNode in self.lastInsertedNodes:
        if lastInsertedNode.insertNode(node):
          nodesWhereInserted.add(lastInsertedNode)
          if node.parent is not None:
            break
      self.lastInsertedNodes = self.lastInsertedNodes - nodesWhereInserted
      self.lastInsertedNodes.add(node)

  def pruneDuplicateNodesOfSameUniProtId(self):
    if self.root:
      for uniprot_id in self.nodes_of_uniprot_id:
        if len(self.nodes_of_uniprot_id[uniprot_id]) > 1:
          all_nodes = list(self.nodes_of_uniprot_id[uniprot_id])
          kept_node = all_nodes[0]
          nodes_to_prune = all_nodes[1:]
          for node in nodes_to_prune:
            if node.parent != kept_node.parent:
              print "Weirdness! Uniprot_id ", uniprot_id, 
              print " is in two different nodes, ", kept_node.left_id,
              print " and ", node.left_id, ", in two different PHOGs, ",
              print kept_node.parent.left_id, " and ", node.parent.left_id,
              print ", in tree ", node.tree_id
            node.parent.children.remove(node)
          self.nodes_of_uniprot_id[uniprot_id] = set([kept_node])

  def computeOrthologs(self):
    if self.root:
      self.root.pruneEmptyChildren()
      self.pruneDuplicateNodesOfSameUniProtId()
      self.root.sortChildren()
      self.uniprot_id_array = numpy.zeros(self.num_uniprot_ids * 2, dtype='i')
      left_right_id_array \
        = numpy.zeros(2 * self.total_num_inserted_nodes, dtype = 'i')
      right_id_of_left_id = {}
      index_of_left_id_in_uniprot_id_array = {}
      node_of_left_id = {}
      node_of_index = {}
      (i_end, k_end) = self.root.fillUniProtIds(self.uniprot_id_array, 
                                self.index_of_uniprot_id, node_of_index, 0,
                                index_of_left_id_in_uniprot_id_array,
                                left_right_id_array, right_id_of_left_id, 
                                node_of_left_id, 0)
      self.end_uniprot_index = i_end
      self.thresholds_of_orthology \
        = numpy.zeros( (i_end, i_end), dtype='d')
      self.phogs_supporting_orthology \
        = numpy.zeros( (i_end, i_end), dtype='i')
      stack_of_containing_phogs = {}
      current_stack_index = -1
      for i in range(i_end):
        query_node = node_of_index[i]
        for k in range(k_end):
          if left_right_id_array[k] in right_id_of_left_id:
            # This is a left_id
            left_id = left_right_id_array[k]
            right_id = right_id_of_left_id[left_id]
            if left_id <= query_node.left_id \
                and right_id >= query_node.right_id:
              node = node_of_left_id[left_id]
              if node.greatest_duplication_distance_of_maximal_descendant:
                # We're entering a PHOG that contains the query sequence, so
                # push it on the stack
                current_stack_index += 1
                stack_of_containing_phogs[current_stack_index] = node
            if left_id in index_of_left_id_in_uniprot_id_array:
              j = index_of_left_id_in_uniprot_id_array[left_id]
              phog = stack_of_containing_phogs[current_stack_index]
              self.phogs_supporting_orthology[i,j] = phog.left_id
              self.thresholds_of_orthology[i,j] = max(0.0, \
                phog.greatest_duplication_distance_of_maximal_descendant)
          else:
            # This is a right_id
            right_id = left_right_id_array[k]
            if current_stack_index >= 0 and right_id == \
                stack_of_containing_phogs[current_stack_index].right_id:
              # We're exiting this PHOG, so pop it off the stack
              del stack_of_containing_phogs[current_stack_index]
              current_stack_index -= 1
      
  def printSelf(self, outf = sys.stdout):
    if self.root:
      self.root.printSelf(level = 0, outf = outf)
      outf.flush()
    else:
      for node in self.lastInsertedNodes:
        root = node.findRoot()
        root.printSelf(level = 0, outf = outf)
      outf.flush()


def serve_trees(comm, tree_server_num, num_tree_servers,
                num_uniprot_processors, all_uniprot_ids):
  served_trees = {}

  tree_row_info = numpy.zeros(6, dtype='d')
  status = MPI.Status()

  while True:
    comm.Recv([tree_row_info, MPI.DOUBLE_PRECISION], source=0, tag=MPI.ANY_TAG,
              status = status)
    if status.Get_tag() == TAG_DATABASE_DONE:
      break
    else:
      tree_id = int(tree_row_info[0])
      try:
        served_tree = served_trees[tree_id]
      except KeyError:
        served_tree = ServedTree()
        served_trees[tree_id] = served_tree
      served_tree.addNodeFromRow(tree_row_info)

  t1 = MPI.Wtime()
  for tree_id in served_trees:
    served_tree = served_trees[tree_id]
    served_tree.computeOrthologs()
    
  t2 = MPI.Wtime()
  print "Tree server %d computed orthologs for %d trees in %g secs" \
    % (tree_server_num, len(served_trees), t2-t1)
  comm.Barrier()

  ortholog_request = numpy.zeros(2, dtype='i')

  while True:
    comm.Recv([ortholog_request, MPI.INT], source=MPI.ANY_SOURCE,
              tag=MPI.ANY_TAG, status = status)
    if status.Get_tag() == TAG_ORTHOLOG_DONE:
      break
    else:
      uniprot_processor_id = status.Get_source()
      tree_id = ortholog_request[0]
      uniprot_id = ortholog_request[1]
      served_tree = served_trees[tree_id]
      if served_tree.root:
        try:
          i = served_tree.index_of_uniprot_id[uniprot_id]
          comm.Send([served_tree.uniprot_id_array[0:
                                            served_tree.end_uniprot_index], 
                    MPI.INT],
                dest=uniprot_processor_id, tag=TAG_ORTHOLOG_RESPONSE)
          comm.Send(
            [served_tree.phogs_supporting_orthology[i,
                                            0:served_tree.end_uniprot_index],
                    MPI.INT],
                  dest=uniprot_processor_id, tag=TAG_ORTHOLOG_RESPONSE)
          comm.Send([served_tree.thresholds_of_orthology[i,
                                            0:served_tree.end_uniprot_index],
                    MPI.DOUBLE_PRECISION],
                  dest=uniprot_processor_id, tag=TAG_ORTHOLOG_RESPONSE)
        except KeyError:
          print "KeyError: tree %d does not have %d" % (tree_id, uniprot_id)
          comm.Send([MPI.BOTTOM, MPI.INT], dest=uniprot_processor_id, tag =
                  TAG_ORTHOLOG_RESPONSE)
          comm.Send([MPI.BOTTOM, MPI.INT], dest=uniprot_processor_id, tag =
                  TAG_ORTHOLOG_RESPONSE)
          comm.Send([MPI.BOTTOM, MPI.DOUBLE_PRECISION], 
                  dest=uniprot_processor_id, tag = TAG_ORTHOLOG_RESPONSE)
      else:
        comm.Send([MPI.BOTTOM, MPI.INT], dest=uniprot_processor_id, tag =
                  TAG_ORTHOLOG_RESPONSE)
        comm.Send([MPI.BOTTOM, MPI.INT], dest=uniprot_processor_id, tag =
                  TAG_ORTHOLOG_RESPONSE)
        comm.Send([MPI.BOTTOM, MPI.DOUBLE_PRECISION], 
                  dest=uniprot_processor_id, tag = TAG_ORTHOLOG_RESPONSE)

def process_uniprot_ids(comm, uniprot_processor_num, num_tree_servers,
                        num_uniprot_processors, all_uniprot_ids):
  trees_of_uniprot_id = {}
  uniprot_ids_in_tree = {}

  uniprot_row_info = numpy.zeros(3, dtype='i')
  status = MPI.Status()

  while True:
    comm.Recv([uniprot_row_info, MPI.INT], source=0, tag=MPI.ANY_TAG,
              status = status)
    if status.Get_tag() == TAG_DATABASE_DONE:
      break
    else:
      tree_id = uniprot_row_info[0]
      left_id = uniprot_row_info[1]
      uniprot_id = uniprot_row_info[2]
      try:
        trees_of_uniprot_id[uniprot_id].add(tree_id)
      except KeyError:
        trees_of_uniprot_id[uniprot_id] = set([tree_id])
      try:
        uniprot_ids_in_tree[tree_id].add(uniprot_id)
      except KeyError:
        uniprot_ids_in_tree[tree_id] = set([uniprot_id])

  comm.Barrier()

  ortholog_request = numpy.zeros(2, dtype='i')
  uniprot_id_array = numpy.zeros(5000, dtype='i')
  phogs_supporting_orthology = numpy.zeros(5000, dtype='i')
  thresholds_of_orthology = numpy.zeros(5000, dtype='d')
  dir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5/'
  f = open(os.path.join(dir, "info_of_uniprot_accession.pkl"))
  info_of_uniprot_accession = cPickle.load(f)
  f.close()
  f = open(os.path.join(dir, "uniprot_accessions_of_uniprot_id.pkl"))
  uniprot_accessions_of_uniprot_id = cPickle.load(f)
  f.close()
  f = open(os.path.join('/clusterfs/vasudha/bpg/OrthologsForQuest/',
                        'OrthologsIn13ReferenceProteomes_%d_of_%d'
                        % (uniprot_processor_num, num_uniprot_processors)), "w")

  base_tree_server_id = 1
  def write_uniprot_id(uniprot_id):
    uniprot_accessions = uniprot_accessions_of_uniprot_id[uniprot_id]
    f.write("%d (%s)" % (uniprot_id, 
                        ','.join(["%s:%s" % (accession,
                        info_of_uniprot_accession[accession]['taxon'])
                        for accession in
                        uniprot_accessions_of_uniprot_id[uniprot_id]])))

  print "UniProt processor %d writing %d uniprot_ids" \
        % (uniprot_processor_num, len(trees_of_uniprot_id))
  t1 = MPI.Wtime()
  for uniprot_id in trees_of_uniprot_id:
    orthologs = {}
    for tree_id in trees_of_uniprot_id[uniprot_id]:
      tree_server_num = base_tree_server_id + tree_id % num_tree_servers
      ortholog_request[0] = tree_id
      ortholog_request[1] = uniprot_id
      comm.Send([ortholog_request, MPI.INT], dest = tree_server_num, 
                tag = TAG_ORTHOLOG_REQUEST)
      comm.Recv([uniprot_id_array, MPI.INT], source = tree_server_num,
                tag = TAG_ORTHOLOG_RESPONSE, status = status)
      num_uniprot_ids = status.Get_count(datatype = MPI.INT)
      comm.Recv([phogs_supporting_orthology, MPI.INT],
                source = tree_server_num, tag = TAG_ORTHOLOG_RESPONSE)
      comm.Recv([thresholds_of_orthology, MPI.DOUBLE_PRECISION], 
                source = tree_server_num, tag = TAG_ORTHOLOG_RESPONSE)
      for i in range(num_uniprot_ids):
        try:
          orthologs[uniprot_id_array[i]].add(
            (tree_id, phogs_supporting_orthology[i],
              thresholds_of_orthology[i]))
        except KeyError:
          orthologs[uniprot_id_array[i]] = set([
            (tree_id, phogs_supporting_orthology[i],
              thresholds_of_orthology[i])])
    write_uniprot_id(uniprot_id)
    f.write(": ")
    for ortholog in orthologs.keys():
      f.write("{")
      write_uniprot_id(ortholog)
      f.write(" <= ")
      for tree_id, left_id, threshold in orthologs[ortholog]:
        f.write("(PHOG%07d_%05d, %f)," % (tree_id, left_id, threshold))
      f.write("};")
    f.write("\n")
  f.close()
  t2 = MPI.Wtime()
  print "UniProt processor %d wrote all uniprot_ids in %g secs" \
        % (uniprot_processor_num, t2-t1)

  comm.Send([MPI.BOTTOM, MPI.INT], dest = 0, tag = TAG_ORTHOLOG_DONE)

def read_from_database(comm, num_tree_servers, num_uniprot_processors, 
                      all_uniprot_ids):
  base_tree_server_id = 1
  base_uniprot_processor_id = 1 + num_tree_servers
  tree_row_info = numpy.zeros(6, dtype='d')
  uniprot_row_info = numpy.zeros(3, dtype='i')

  bpg_password = get_credentials('bpg_user')
  connection = psycopg2.connect(
    "dbname='%s' user='%s' host='db' password='%s'" %
    ('pfacts003_test', 'bpg_user', bpg_password))
  cur = connection.cursor('ortholog_cursor',
                          cursor_factory = psycopg2.extras.DictCursor)

  db_row_fetch_time = 0.0
  tree_row_send_time = 0.0
  uniprot_row_send_time = 0.0
  row_prep_time = 0.0
  t1 = MPI.Wtime()
  sql = """SELECT tree_id,
                  tree_node_left_id,
                  tree_node_right_id,
                  duplication_distance,
                  greatest_duplication_distance_of_maximal_descendant,
                  uniprot_id
                  FROM tree_node_uniprot_taxonomy_materialized
                  """

  cur.execute(sql)
  num_database_rows = 0
  db_row_start_t = MPI.Wtime()
  for row in cur:
    db_row_end_t = MPI.Wtime()
    db_row_fetch_time += db_row_end_t - db_row_start_t
    cur_tree_row_send_time = 0.0
    cur_uniprot_row_send_time = 0.0
    num_database_rows += 1
    if num_database_rows > 0 and num_database_rows % 1000000 == 0:
      print "Read %d rows from database so far" % num_database_rows
    tree_id = row[0]
    greatest_duplication_distance_of_maximal_descendant = row[4]
    uniprot_id = row[5]
    if greatest_duplication_distance_of_maximal_descendant and \
          (not row[3] or row[3] > row[4]) or \
        uniprot_id is not None and uniprot_id in all_uniprot_ids:
      for i in range(6):
        if row[i] is not None:
          tree_row_info[i] = float(row[i])
        else:
          tree_row_info[i] = None
      if uniprot_id is not None and uniprot_id not in all_uniprot_ids:
        tree_row_info[5] = None
      tree_server_num = tree_id % num_tree_servers
      tree_row_send_start_t = MPI.Wtime()
      comm.Send([tree_row_info,MPI.DOUBLE_PRECISION], 
                dest=tree_server_num + base_tree_server_id,
                tag=TAG_DATABASE_ROW)
      tree_row_send_end_t = MPI.Wtime()
      cur_tree_row_send_time = tree_row_send_end_t - tree_row_send_start_t
      tree_row_send_time += cur_tree_row_send_time
      if not greatest_duplication_distance_of_maximal_descendant:
        uniprot_row_info[0] = int(row[0]) # tree_id
        uniprot_row_info[1] = int(row[1]) # tree_node_left_id
        uniprot_row_info[2] = int(row[5]) # uniprot_id
        uniprot_processor_num = uniprot_id % num_uniprot_processors
        uniprot_row_send_start_t = MPI.Wtime()
        comm.Send([uniprot_row_info,MPI.INT],
                  dest=uniprot_processor_num + base_uniprot_processor_id,
                  tag=TAG_DATABASE_ROW)
        uniprot_row_send_end_t = MPI.Wtime()
        cur_uniprot_row_send_time = uniprot_row_send_end_t \
                                    - uniprot_row_send_start_t
        uniprot_row_send_time += cur_uniprot_row_send_time
    db_row_start_t = MPI.Wtime()
    row_prep_time += (db_row_start_t - db_row_end_t) \
                  - cur_tree_row_send_time - cur_uniprot_row_send_time

  t2 = MPI.Wtime()
  print "Finished reading ", num_database_rows, 
  print " rows of the database in ", t2 - t1, " secs"
  print "Total time fetching database rows: ", db_row_fetch_time
  print "Total time sending tree rows: ", tree_row_send_time
  print "Total time sending uniprot rows: ", uniprot_row_send_time
  print "Total time preparing rows to send out: ", row_prep_time
  for tree_server_num in range(num_tree_servers):
    comm.Send([MPI.BOTTOM,MPI.INT], dest=tree_server_num + base_tree_server_id, 
              tag=TAG_DATABASE_DONE)

  for uniprot_processor_num in range(num_uniprot_processors):
    comm.Send([MPI.BOTTOM,MPI.INT], 
              dest=uniprot_processor_num + base_uniprot_processor_id,
              tag=TAG_DATABASE_DONE)
    
def coordinate_termination(comm, num_tree_servers, num_uniprot_processors):
  t1 = MPI.Wtime()
  base_tree_server_id = 1
  status = MPI.Status()
  for i in range(num_uniprot_processors):
    comm.Recv([MPI.BOTTOM,MPI.INT], source = MPI.ANY_SOURCE, tag =
                TAG_ORTHOLOG_DONE, status = status)
  for tree_server_num in range(num_tree_servers):
    comm.Send([MPI.BOTTOM,MPI.INT], 
              dest = tree_server_num + base_tree_server_id,
              tag = TAG_ORTHOLOG_DONE)
  t2 = MPI.Wtime()
  print "Finished reading and writing orthologs in ", t2 - t1, " secs"

def main():
  comm = MPI.COMM_WORLD
  myid = comm.Get_rank()
  numprocs = comm.Get_size()

  num_tree_servers = (numprocs - 1)/ 2
  num_uniprot_processors = (numprocs - 1)- num_tree_servers
  linenum = 0
  dir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5'
  os.chdir(dir)
  all_uniprot_ids = set()
  f = open('all_uniprot_ids.txt')
  for line in f.readlines():
    all_uniprot_ids.add(int(line.strip()))
  f.close()
  if myid == 0:
    read_from_database(comm, num_tree_servers, num_uniprot_processors,
                      all_uniprot_ids)
    t1 = MPI.Wtime()
    comm.Barrier()
    t2 = MPI.Wtime()
    print "All tree servers have finished computing orthologs in %g secs" \
          % (t2 - t1)
    coordinate_termination(comm, num_tree_servers, num_uniprot_processors)
  elif myid - 1 < (numprocs - 1)/ 2:
    serve_trees(comm, myid - 1, num_tree_servers, num_uniprot_processors,
                all_uniprot_ids)
  else:
    process_uniprot_ids(comm, myid - 1 -  (numprocs - 1)/ 2, num_tree_servers, 
                        num_uniprot_processors, all_uniprot_ids)

if __name__ == '__main__':
  main()
