#!/usr/bin/env python

import os
from pfacts003.phylofacts.models import Tree, TreeNode
from dir_of_family import get_dir_of_family_accession
from array import array
from optparse import OptionParser

def tree_escape(buffer):
  return buffer.replace('(', '%28').replace(')', '%29').replace(':',
                              '%3a').replace(',', '%2c').replace(';', '%3b')

def get_summary_array(summary_array, nodes, nodeIndex, threshold, printing):
  node = nodes[nodeIndex]
  if node.sequence_header is not None:
    if printing:
      summary_array.fromstring("SEQHDR%d" % node.sequence_header.id)
      if node.distance_to_parent is not None and node.distance_to_parent >= 0:
        summary_array.fromstring(':%g' % node.distance_to_parent)
    printing = False
  elif node.greatest_duplication_distance_of_maximal_descendant is not None \
      and node.greatest_duplication_distance_of_maximal_descendant < threshold \
      and (node.duplication_distance is None 
            or node.duplication_distance >= threshold):
    if printing:
      summary_array.fromstring("PHOG%07d_%05d" % (node.tree.id, node.left_id))
      if node.distance_to_parent is not None and node.distance_to_parent >= 0:
        summary_array.fromstring(':%g' % node.distance_to_parent)
    printing = False
  if printing:
    # Add the '(' to the string
    summary_array.fromstring('(')
  nodeIndex += 1
  firstChild = True
  while nodeIndex < len(nodes) and nodes[nodeIndex].right_id < node.right_id:
    if not firstChild and printing:
      summary_array.fromstring(',')
    nodeIndex = get_summary_array(summary_array, nodes, nodeIndex,
                                      threshold, printing)
    firstChild = False
  if printing:
    summary_array.fromstring(')')
    if node.likelihood is not None:
      summary_array.fromstring('%g' % node.likelihood)
    if node.distance_to_parent is not None and node.distance_to_parent >= 0:
      summary_array.fromstring(':%g' % node.distance_to_parent)
  return nodeIndex

def phog_summary_tree(family_accession, method='ml', threshold=0.0):
  tree = Tree.objects.get(family__id = int(family_accession[3:]), method=method)
  nodes = TreeNode.objects.filter(tree = tree).order_by('left_id')
  summary_array = array('c')
  get_summary_array(summary_array, nodes, 0, threshold, True)
  summary_array.fromstring(';')
  return summary_array.tostring()

def main():
  usage = "%prog [options] family_accession"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-t", "--threshold", dest="threshold", default=0.0,
                    help="PHOG threshold at which to cut tree")
  opt_parser.add_option("-m", "--method", dest="method", default='ml',
                    help="Tree method")
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  try:
    threshold = float(options.threshold)
  except ValueError:
    opt_parser.error('Threshold must be a nonnegative real number.')
  if threshold < 0.0:
    opt_parser.error('Threshold must be a nonnegative real number.')
  family_accession = args[0]
  summary_tree = phog_summary_tree(family_accession, options.method, threshold)
  outpath = os.path.join(get_dir_of_family_accession(family_accession), 
            '%s_phogt_%g.%s' % (family_accession, threshold, options.method))
  outf = open(outpath, 'w')
  outf.write(summary_tree)
  outf.write('\n')
  outf.close()

if __name__ == '__main__':
  main()
