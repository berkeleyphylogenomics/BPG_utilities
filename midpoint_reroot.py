#!/usr/bin/env python2.6

from cogent.core.tree import PhyloNode
from cogent import LoadTree
from cogent.parse.tree import DndParser
from optparse import OptionParser
import os, sys
import cProfile
from datetime import date
import time
import numpy

def get_breadth_first_visit_order(tree):
  breadth_first_visit_order = {}
  visit_order_of_node = {}
  current_index = 0
  num_nodes = 0
  breadth_first_visit_order[0] = tree
  visit_order_of_node[tree] = 0
  num_nodes += 1
  for child in tree.Children:
    breadth_first_visit_order[num_nodes] = child
    visit_order_of_node[child] = num_nodes
    num_nodes += 1
  current_index += 1
  while current_index < num_nodes:
    current_node = breadth_first_visit_order[current_index]
    for child in current_node.Children:
      breadth_first_visit_order[num_nodes] = child
      visit_order_of_node[child] = num_nodes
      num_nodes += 1
    current_index += 1
  child_visit_orders_of = {}
  branch_length_of = numpy.zeros(num_nodes)
  for i in xrange(num_nodes):
    node = breadth_first_visit_order[i]
    branch_length_of[i] = node.params['length']
    child_visit_orders_of[i] = set()
    for child in node.Children:
      child_visit_orders_of[i].add(visit_order_of_node[child])
  return (breadth_first_visit_order, visit_order_of_node, 
          branch_length_of, child_visit_orders_of, num_nodes)

def getTimeStr(starttime, endtime):

    """Returns time interval formatted in minutes and seconds

    This function formats a time interval into minutes and seconds, so we can
display the elapsed time for users
"""
    elapsed_secs = endtime - starttime

    if elapsed_secs > 60:
        mins = int(elapsed_secs / 60)
        elapsed_secs = elapsed_secs - (mins * 60)

        if mins == 1:
            return '[1 min, %.2f secs]' % elapsed_secs

        return '[%d mins, %.2f secs]' % (mins, elapsed_secs)

    return '[%.2f secs]' % elapsed_secs


def main():
  usage = "%prog [options] tree_to_midpoint_reroot"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()

  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if not os.path.exists(args[0]):
    opt_parser.error('Tree file %s not found' % args[0])

  f = open(args[0])
  tree_string = f.read()
  f.close()

  unrooted_tree = DndParser(tree_string, PhyloNode)
  breadth_first_visit_order, visit_order_of_node, branch_length_of, \
      child_visit_orders_of, num_nodes \
      = get_breadth_first_visit_order(unrooted_tree)

  # We will refer to the node objects by their index in the visit order

  tip_node_objects = unrooted_tree.tips()
  num_tips = len(tip_node_objects)

  # We will refer to the tip objects by their index in the tip_node_objects
  # list

  distance_from_node_to_tip = numpy.zeros((num_nodes, num_tips))
  stepping_stone_from_node_to_tip = numpy.zeros((num_nodes, num_tips))
  tips_connected_to_node = {}
  
  for node in xrange(num_nodes):
    tips_connected_to_node[node] = set()
    for tip in xrange(num_tips):
      distance_from_node_to_tip[node,tip] = -1.0
      stepping_stone_from_node_to_tip[node,tip] = -1

  for tip in xrange(num_tips):
    tip_as_node = visit_order_of_node[tip_node_objects[tip]]
    distance_from_node_to_tip[tip_as_node,tip] = 0.0
    stepping_stone_from_node_to_tip[tip_as_node,tip] = tip_as_node
    tips_connected_to_node[tip_as_node].add(tip)

  for parent in reversed(xrange(num_nodes)):
    for child in child_visit_orders_of[parent]:
      child_to_parent_distance = branch_length_of[child]
      for tip in tips_connected_to_node[child]:
        tip_distance_to_child = distance_from_node_to_tip[child, tip]
        tip_to_parent_distance_through_child \
          = tip_distance_to_child + child_to_parent_distance
        if tip in tips_connected_to_node[parent]:
          tip_distance_to_parent = distance_from_node_to_tip[parent, tip]
          if tip_to_parent_distance_through_child < tip_distance_to_parent:
            distance_from_node_to_tip[parent, tip] \
                = tip_to_parent_distance_through_child
            stepping_stone_from_node_to_tip[parent, tip] = child
        else:
          distance_from_node_to_tip[parent, tip] \
              = tip_to_parent_distance_through_child
          stepping_stone_from_node_to_tip[parent, tip] = child
          tips_connected_to_node[parent].add(tip)

  for parent in xrange(num_nodes):
    for child in child_visit_orders_of[parent]:
      child_to_parent_distance = branch_length_of[child]
      for tip in tips_connected_to_node[parent]:
        tip_distance_to_parent = distance_from_node_to_tip[parent, tip]
        tip_to_child_distance_through_parent \
          = tip_distance_to_parent + child_to_parent_distance
        if tip in tips_connected_to_node[child]:
          tip_distance_to_child = distance_from_node_to_tip[child, tip]
          if tip_to_child_distance_through_parent < tip_distance_to_child:
            distance_from_node_to_tip[child, tip] \
                = tip_to_child_distance_through_parent
            stepping_stone_from_node_to_tip[child, tip] = parent
        else:
          distance_from_node_to_tip[child, tip] \
              = tip_to_child_distance_through_parent
          stepping_stone_from_node_to_tip[child, tip] = parent
          tips_connected_to_node[child].add(tip)

  max_distance = 0.0
  max_tip0 = None
  max_tip1 = None

  for tip0 in xrange(num_tips):
    tip0_as_node = visit_order_of_node[tip_node_objects[tip0]]
    for tip1 in xrange(num_tips):
      tip0_tip1_distance = distance_from_node_to_tip[tip0_as_node,tip1]
      if tip0_tip1_distance > max_distance:
        max_distance = tip0_tip1_distance
        max_tip0 = tip0
        max_tip1 = tip1

  midpoint_distance = max_distance / 2
  tip0 = max_tip0
  tip1 = max_tip1
  node_closer_to_tip0 = visit_order_of_node[tip_node_objects[tip0]]
  node_closer_to_tip1 = node_closer_to_tip0
  distance_to_tip1 = distance_from_node_to_tip[node_closer_to_tip0, tip1]
  node_even_closer_to_tip1 \
      = stepping_stone_from_node_to_tip[node_closer_to_tip0, tip1]
  previous_distance_to_tip1 = distance_to_tip1
  while distance_to_tip1 > midpoint_distance:
    node_closer_to_tip0 = node_closer_to_tip1
    node_closer_to_tip1 = node_even_closer_to_tip1
    previous_distance_to_tip1 = distance_to_tip1
    distance_to_tip1 = distance_from_node_to_tip[node_closer_to_tip1, tip1]
    node_even_closer_to_tip1 \
        = stepping_stone_from_node_to_tip[node_closer_to_tip1, tip1]

  node_object_closer_to_tip0 = breadth_first_visit_order[node_closer_to_tip0]
  node_object_closer_to_tip1 = breadth_first_visit_order[node_closer_to_tip1]
  if node_object_closer_to_tip1 == node_object_closer_to_tip0._parent:
    theParent = node_object_closer_to_tip1
    theChild = node_object_closer_to_tip0
    distance_from_new_root_to_parent \
      = midpoint_distance - distance_to_tip1
    distance_from_new_root_to_child \
      = previous_distance_to_tip1 - midpoint_distance
  elif node_object_closer_to_tip0 == node_object_closer_to_tip1._parent:
    theParent = node_object_closer_to_tip0
    theChild = node_object_closer_to_tip1
    distance_from_new_root_to_parent \
      = previous_distance_to_tip1 - midpoint_distance
    distance_from_new_root_to_child \
      = midpoint_distance - distance_to_tip1
  else:
    # Should never get here
    raise AssertionError('Adjacent nodes on maximum span not parent-child')

  sys.stdout.write('(')
  # omit the branch length from theChild to its parent, since this is the
  # branch being broken in two
  sys.stdout.write(':'.join(
    theChild.getNewick(with_distances=True,semicolon=False).split(':')[:-1]))
  sys.stdout.write(":%g" % distance_from_new_root_to_child)
  sys.stdout.write(',')
  def print_rotated_node(child, parent, distance_from_parent_to_new_parent):
    sys.stdout.write('(')
    sys.stdout.write(','.join([other_child.getNewick(with_distances=True,
                                                    semicolon=False)
                                for other_child in parent.Children
                                if other_child != child]))
    if parent._parent:
      sys.stdout.write(',')
      print_rotated_node(parent, parent._parent, parent.params['length'])
    sys.stdout.write(')')
    if parent.Name:
      sys.stdout.write(parent.Name)
    sys.stdout.write(":%g" % distance_from_parent_to_new_parent)
  print_rotated_node(theChild, theParent, distance_from_new_root_to_parent)
  sys.stdout.write(');\n')

if __name__ == '__main__':
  main()
