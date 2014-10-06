#!/usr/bin/env python

import sys
from pfacts003.phylofacts.models import Tree, TreeNode, SequenceHeader, \
    TreeNodeAlignment, OrthologTypes, preset_thresholds, \
    get_ortholog_type_of_threshold
from optparse import OptionParser

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

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession sequence_identifier"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-m", "--method", dest="method", default="ml",
                      help="tree method for which to extract orthologs.")
  opt_parser.add_option("-t", "--threshold", dest="threshold",
              default="0.0",
              help="Threshold at which to extract orthologs.")
  opt_parser.add_option("-d", "--include_descriptions",
              action="store_true", default=False,
              help="Whether to include the full descriptions in the tree")
  opt_parser.add_option("-v", "--verbose", dest="verbose",
              action="store_true", default=True,
              help="Whether to print verbose output")
  opt_parser.add_option("-q", "--quiet", dest="verbose",
              action="store_false", default=True,
              help="Whether to suppress verbose output")
  (options, args) = opt_parser.parse_args()
  if len(args) < 2:
    opt_parser.error("Incorrect number of arguments")
  else:
    bpg_accession = args[0]
    sequence_header_identifier = args[1].split()[0]
  try:
    opt_threshold = float(options.threshold)
  except ValueError:
    opt_parser.error("Threshold must be a non-negative floating-point number")
  if opt_threshold < 0.0:
    opt_parser.error("Threshold must be a non-negative floating-point number")
  verbose = options.verbose
  tree = get_tree(bpg_accession, options.method)
  if not tree:
    opt_parser.error('Unable to find %s tree for book %s in database'
                      % (options.method, bpg_accession))
  leaves = TreeNode.objects.filter(tree = tree, 
                            sequence_header__isnull = False)
  query_leaf = None
  for leaf in leaves:
    if leaf.sequence_header.header.split()[0] == sequence_header_identifier:
      query_leaf = leaf
      break

  if not query_leaf:
    opt_parser.error('Unable to find %s in tree' % sequence_header_identifier)

  threshold = opt_threshold
  ortholog_type = get_ortholog_type_of_threshold(threshold)

  phog = query_leaf.get_containing_phog(threshold=threshold)

  if not phog:
    print "The sequence is not in a PHOG at this threshold."
    sys.exit(0)

  phog_accession \
    = phog.get_accession(ortholog_type = ortholog_type, threshold=threshold)
  if verbose:
    print phog_accession
  orthologs = phog.get_included_leaves(threshold=threshold).order_by('left_id')
  ortholog_of_sequence_header = {}
  for ortholog in orthologs:
    print ortholog.sequence_header.header
    ortholog_of_sequence_header[ortholog.sequence_header] = ortholog
    
  if verbose:
    print "Writing alignment to %s.a2m..." % phog_accession,
    sys.stdout.flush()
  root_node = TreeNode.objects.filter(tree = query_leaf.tree, left_id = 1)
  tree_node_alignment_objects = TreeNodeAlignment.objects.filter(
    tree_node = root_node,
    sequence_header__in = [ortholog.sequence_header 
                            for ortholog in orthologs])
  alignment_of_ortholog = {}
  for obj in tree_node_alignment_objects:
    alignment_of_ortholog[ortholog_of_sequence_header[obj.sequence_header]] \
        = obj

  f = open('%s.a2m' % phog_accession,"w")
  for ortholog in orthologs:
    if ortholog not in alignment_of_ortholog:
      continue
    obj = alignment_of_ortholog[ortholog]
    f.write('>%s\n' % obj.sequence_header.header)
    f.write('%s\n' % obj.aligned_sequence.chars)
  f.close()
  if verbose:
    print ("done.")
    sys.stdout.flush()

  if verbose:
    print "Writing tree to %s.newick..." % phog_accession,
    sys.stdout.flush()
  included_tree_nodes = set()
  for ortholog in orthologs:
    intermediate_tree_nodes = TreeNode.objects.filter(tree = tree,
                                left_id__lte = ortholog.left_id,
                                right_id__gte = ortholog.right_id,
                                left_id__gte = phog.left_id,
                                right_id__lte = phog.right_id)
    for node in intermediate_tree_nodes:
      included_tree_nodes.add(node)
  left_ids = set([node.left_id for node in included_tree_nodes])
  children_of_left_id = {}
  for left_id in left_ids:
    children_of_left_id[left_id] = set()
  for node in included_tree_nodes:
    if node.parent_left_id in left_ids:
      children_of_left_id[node.parent_left_id].add(node)
  f = open('%s.newick' % phog_accession,"w")
  def write_tree_node(node):
    if node.sequence_header:
      if options.include_descriptions:
        escaped_header = node.sequence_header.header
        escaped_header = escaped_header.replace('(', '%28')
        escaped_header = escaped_header.replace(')', '%29')
        escaped_header = escaped_header.replace(':', '%3A')
        escaped_header = escaped_header.replace(',', '%2C')
        escaped_header = escaped_header.replace(';', '%3B')
        f.write(escaped_header)
      else:
        f.write('%s' % node.sequence_header.header.split()[0])
    else:
      f.write('(')
      first_child = True
      for child in children_of_left_id[node.left_id]:
        if first_child:
          first_child = False
        else:
          f.write(',')
        write_tree_node(child)
      f.write(')')
    if node.likelihood:
      f.write('%f' % node.likelihood)
    if node.distance_to_parent:
      f.write(':%f' % node.distance_to_parent)
  write_tree_node(phog)
  f.write(';\n')
  f.close()
  if verbose:
    print "done."
    sys.stdout.flush()

if __name__ == '__main__':
  main()
