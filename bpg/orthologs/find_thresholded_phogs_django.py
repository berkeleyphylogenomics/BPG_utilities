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

from optparse import OptionParser
from django.db.models import Q
from pfacts003.phylofacts.models import Tree, SpeciesNearestMaximalTreeNode, \
  TreeNode
from django.db import transaction

def get_duplication_distance(tree_node_id):
  maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                                          tree_node__id__exact = tree_node_id)
  if maximal_tree_nodes:
    maximal_tree_node = maximal_tree_nodes[0]
    if maximal_tree_node.next_nearest_distance >= 0.0:
      return (maximal_tree_node.distance + \
              maximal_tree_node.next_nearest_distance) / 2
    else:
      return 5000
  else:
    return -1.0

def get_difference_distance(tree_node_id):
  maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                                          tree_node__id__exact = tree_node_id)
  if maximal_tree_nodes:
    maximal_tree_node = maximal_tree_nodes[0]
    if maximal_tree_node.next_nearest_distance >= 0.0:
      return (maximal_tree_node.next_nearest_distance - \
              maximal_tree_node.next_nearest_distance)
    else:
      return 5000
  else:
    return -1.0

class ProximalTreeNode:
  def __init__(self, visit_id, tree_node_id, left_id, right_id, 
                parent_visit_id, verbose=False):
    if visit_id > 0:
      duplication_distance = get_duplication_distance(tree_node_id)
    else:
      duplication_distance = 5000.0
    if visit_id > 0:
      difference_distance = get_difference_distance(tree_node_id)
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

def get_tree(bpg_accession, method='nj'):
  try:
    family_id = int(bpg_accession[3:])
  except ValueError:
    return None
  tree_objects = Tree.objects.filter(family__id__exact = family_id,
                                      method__exact = method)
  if tree_objects:
    return tree_objects[0]
  else:
    return None

def print_all_proximal_nodes(tree):
  maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                            tree_node__tree__exact =
                            tree).order_by('tree_node.left_id')
  print "%d proximal nodes" % maximal_tree_nodes.count()
  print "-------------------------------------"
  print "| tree_node_id | left_id | right_id |"
  print "-------------------------------------"
  for maximal_tree_node in maximal_tree_nodes:
    tree_node = maximal_tree_node.tree_node
    print "|    %09d | %07d |  %07d |" \
      % (tree_node.id, tree_node.left_id, tree_node.right_id)
    print "-------------------------------------"
  print
    
def get_breadth_first_proximal_nodes(tree, verbose=False):
  breadth_first_search_order = {}
  current_index = 0
  num_nodes = 0
  fake_root = ProximalTreeNode(num_nodes, -1, -1, -1, -1, verbose)
  breadth_first_search_order[0] = fake_root
  num_nodes += 1
  maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                            tree_node__tree__exact =
                            tree).order_by('tree_node.left_id')[0:1]
  while maximal_tree_nodes:
    maximal_tree_node = maximal_tree_nodes[0]
    tree_node = maximal_tree_node.tree_node
    breadth_first_search_order[num_nodes] \
      = ProximalTreeNode(num_nodes, tree_node.id, 
                          tree_node.left_id, tree_node.right_id, 0, verbose)
    fake_root.children_visit_ids.add(num_nodes)
    num_nodes += 1
    maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                        tree_node__tree__exact = tree,
                        tree_node__left_id__gt =
                          tree_node.right_id).order_by(
                                                    'tree_node.left_id')[0:1]
  current_index += 1
  while current_index < num_nodes:
    current_node = breadth_first_search_order[current_index]
    maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                        tree_node__tree__exact = tree,
                        tree_node__left_id__gt = current_node.left_id,
                        tree_node__right_id__lt = 
                        current_node.right_id).order_by(
                                                    'tree_node.left_id')[0:1]
    while maximal_tree_nodes:
      maximal_tree_node = maximal_tree_nodes[0]
      tree_node = maximal_tree_node.tree_node
      breadth_first_search_order[num_nodes] \
        = ProximalTreeNode(num_nodes, tree_node.id,
                            tree_node.left_id, tree_node.right_id,
                            current_index, verbose)
      current_node.children_visit_ids.add(num_nodes)
      num_nodes += 1
      maximal_tree_nodes = SpeciesNearestMaximalTreeNode.objects.filter(
                        tree_node__tree__exact = tree,
                        tree_node__left_id__gt = tree_node.right_id,
                        tree_node__right_id__lt = 
                        current_node.right_id).order_by(
                                                'tree_node.left_id')[0:1]
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
              
def make_distances_monotonic(breadth_first_search_order, num_nodes,
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

@transaction.commit_on_success
def put_distances_in_tree_node_table(breadth_first_search_order,
                                    num_nodes, verbose = False):
  # Skip the fake root node at index 0
  for i in xrange(1,num_nodes):
    node = breadth_first_search_order[i]
    if verbose:
      print "Updating duplication distance for node %d to %g" \
          % (node.tree_node_id, node.duplication_distance)
      print "Updating difference distance for node %d to %g" \
          % (node.tree_node_id, node.difference_distance)
    tree_nodes = TreeNode.objects.filter(id = node.tree_node_id)
    tree_nodes.update(duplication_distance = node.duplication_distance,
                      difference_distance = node.difference_distance)

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

@transaction.commit_on_success
def put_phogts_in_tree_node_table(tree, breadth_first_search_order, 
                          num_nodes, phogts, threshold, kind, verbose = False,
                          set_is_phogt_false = False):
  root_tree_node = TreeNode.objects.get(tree__exact = tree, 
                                        left_id__exact = 1)
  if verbose:
    print "Root node id: ", root_tree_node.id
    print "Set phogt_%s_node to %d for all leaves in tree %d" \
        % (kind, root_tree_node.id, tree.id)
  phogt_of_leaf = {}
  leaves = TreeNode.objects.filter(tree__exact = tree, 
                                sequence_header__isnull = False)
  for leaf in leaves:
    phogt_of_leaf[leaf] = root_tree_node.id
  for i in xrange(1, num_nodes):
    if i in phogts:
      node = breadth_first_search_order[i]
      if verbose:
        print "Set phogt_%s_node to %d for all leaves in tree %d " \
          % (kind, node.tree_node_id, tree.id),
        print " with left_id >= %d and right_id <= %d " \
          % (node.left_id, node.right_id)
      for leaf in leaves:
        if leaf.left_id >= node.left_id and leaf.right_id <= node.right_id:
          phogt_of_leaf[leaf] = node.tree_node_id
  leaves_of_phogt = {}
  for leaf in phogt_of_leaf.keys():
    if phogt_of_leaf[leaf] not in leaves_of_phogt:
      leaves_of_phogt[phogt_of_leaf[leaf]] = set()
    leaves_of_phogt[phogt_of_leaf[leaf]].add(leaf.id)
  for phogt in leaves_of_phogt.keys():
    phogt_node = TreeNode.objects.get(id = phogt)
    its_leaves = TreeNode.objects.filter(id__in = leaves_of_phogt[phogt])
    if kind == 'tight':
      its_leaves.update(phogt_tight_node = phogt_node)
    elif kind == 'medium':
      its_leaves.update(phogt_medium_node = phogt_node)
    elif kind == 'loose':
      its_leaves.update(phogt_loose_node = phogt_node)
  tree_node_ids = [breadth_first_search_order[i].tree_node_id for i in phogts]
  if verbose:
    print "Set is_phogt_%s to true for %s" \
        % (kind, ','.join([str(id) for id in tree_node_ids]))
  phogt_nodes = TreeNode.objects.filter(id__in = tree_node_ids)
  if kind == 'tight':
    phogt_nodes.update(is_phogt_tight = True)
  elif kind == 'medium':
    phogt_nodes.update(is_phogt_medium = True)
  elif kind == 'loose':
    phogt_nodes.update(is_phogt_loose = True)
  if set_is_phogt_false:
    if verbose:
      print "Set is_phogt_%s to false for all other nodes in tree %d" \
          % (kind, tree.id)
    all_other_nodes = TreeNode.objects.filter(
                            tree__exact = tree).exclude(id__in = tree_node_ids)
    if kind == 'tight':
      all_other_nodes.update(is_phogt_tight = False)
    elif kind == 'medium':
      all_other_nodes.update(is_phogt_medium = False)
    elif kind == 'loose':
      all_other_nodes.update(is_phogt_loose = False)

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
  tree = get_tree(bpg_accession, options.method)
  if not tree:
    opt_parser.error('Unable to find %s tree for book %s in database'
                      % (options.method, bpg_accession))
  (breadth_first_search_order, num_proximal_nodes) \
    = get_breadth_first_proximal_nodes(tree, verbose)
  make_distances_monotonic(breadth_first_search_order, num_proximal_nodes,
                            verbose)
  put_distances_in_tree_node_table(breadth_first_search_order, 
                            num_proximal_nodes, verbose)
  thresholds = {}
  thresholds['tight'] = 0.09375
  thresholds['medium'] = 0.296874
  thresholds['loose'] = 0.9375
  for kind in thresholds.keys():
    phogts = compute_phogts_at_threshold(breadth_first_search_order,
                          num_proximal_nodes, thresholds[kind], kind, verbose)
    put_phogts_in_tree_node_table(tree, breadth_first_search_order, 
                                  num_proximal_nodes, phogts, thresholds[kind],
                                  kind, verbose)

if __name__ == '__main__':
  main()
