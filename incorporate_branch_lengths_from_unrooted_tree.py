#!/usr/bin/env python

import sys, string, re

id_re = re.compile('SEQ[0-9]*')
branch_length_re = re.compile('[0-9]\.[0-9e\-]*')

class node:
  def __init__(self, seqid=None):
    self.parent = None
    self.children = set()
    self.branch_length = -1.0
    self.seqid = seqid
    if seqid:
      self.contained_seqids = set([seqid])
    else:
      self.contained_seqids = set()

  def addChild(self, child):
    self.children.add(child)
    child.parent = self
    self.contained_seqids = self.contained_seqids | child.contained_seqids

  def readFromTreeString(self, tree_string, i):
    if i >= len(tree_string):
      return len(tree_string)
    if tree_string[i] == '(':
      while tree_string[i] != ')':
        child = node()
        i = child.readFromTreeString(tree_string, i+1)
        self.addChild(child)
      i += 1
    else:
      m = id_re.match(tree_string[i:])
      seqid = m.group(0)
      i += len(seqid)
      self.seqid = seqid
      self.contained_seqids = set([self.seqid])
    if i < len(tree_string) and tree_string[i] == ':':
      i += 1
      branch_length = branch_length_re.match(tree_string[i:]).group(0)
      self.branch_length = float(branch_length)
      i += len(branch_length)
    return i

  def printSelf(self, indentLevel = 0):
    if self.seqid:
      sys.stdout.write("%s" % self.seqid)
    else:
      sys.stdout.write("(")
      children = list(self.children)
      for i in range(len(children)):
        children[i].printSelf(indentLevel + 1)
        if i < len(children) - 1:
          sys.stdout.write(",")
      sys.stdout.write(")")
    if self.branch_length >= 0.0:
      sys.stdout.write(":%g" % self.branch_length)

  def getContainedLeaves(self):
    if self.seqid:
      return set([self])
    else:
      contained_leaves = set()
      for child in self.children:
        contained_leaves = contained_leaves | child.getContainedLeaves()
      return contained_leaves

def get_breadth_first_search_order(root):
  breadth_first_search_order = {}
  breadth_first_search_order[0] = root
  num_nodes = 1
  node_index = 0
  while node_index < num_nodes:
    node = breadth_first_search_order[node_index]
    for child in node.children:
      breadth_first_search_order[num_nodes] = child
      num_nodes += 1
    node_index += 1
  return (num_nodes, breadth_first_search_order)

def incorporate_branch_lengths_from_unrooted_tree(satchmo_root, raxml_root):
  all_seqids = satchmo_root.contained_seqids
  if all_seqids != raxml_root.contained_seqids:
    print "These trees are not the same! Exiting..."
    return
  (num_satchmo_nodes, breadth_first_satchmo_search_order) = \
      get_breadth_first_search_order(satchmo_root)
  raxml_leaves = raxml_root.getContainedLeaves()
  raxml_leaf_of_seqid = {}
  for leaf in raxml_leaves:
    raxml_leaf_of_seqid[leaf.seqid] = leaf
  raxml_branch_of_satchmo_branch = {}
  satchmo_branch_of_raxml_branch = {}
  def match_satchmo_branch_with_raxml_branch(
      satchmo_node, raxml_node):
    raxml_branch_of_satchmo_branch[satchmo_node] = raxml_node
    satchmo_branch_of_raxml_branch[raxml_node] = satchmo_node
    satchmo_node.branch_length = raxml_node.branch_length
  reverse_range = range(1,num_satchmo_nodes)
  reverse_range.reverse()
  for i in reverse_range:
    satchmo_node = breadth_first_satchmo_search_order[i]
    if satchmo_node.seqid:
      raxml_node = raxml_leaf_of_seqid[satchmo_node.seqid]
      match_satchmo_branch_with_raxml_branch(satchmo_node, raxml_node)
    else:
      for satchmo_child in satchmo_node.children:
        if satchmo_child in raxml_branch_of_satchmo_branch:
          raxml_child = raxml_branch_of_satchmo_branch[satchmo_child]
          if raxml_child.parent:
            raxml_node = raxml_child.parent
            if raxml_node.contained_seqids == satchmo_node.contained_seqids:
              match_satchmo_branch_with_raxml_branch(satchmo_node, raxml_node)
              break
  # At this point, all the SATCHMO subtrees which correspond to RAxML subtrees
  # have branch lengths.  Now for the rest of the SATCHMO branches, the
  # corresponding RAxML branches may lie in any direction (i.e., they may be
  # going in the reverse direction from the SATCHMO branches).  So now we have
  # to look in all directions (at both parents and children) to find a node
  # with a matching split.  The split may present as its complement, and it may
  # be on the other end of the branch.
  for i in reverse_range:
    satchmo_node = breadth_first_satchmo_search_order[i]
    if satchmo_node not in raxml_branch_of_satchmo_branch:
      if satchmo_node.parent and \
          satchmo_node.parent in raxml_branch_of_satchmo_branch:
        satchmo_parent = satchmo_node.parent
        raxml_parent = raxml_branch_of_satchmo_branch[satchmo_parent]
        if raxml_parent.parent \
            and raxml_parent.parent not in satchmo_branch_of_raxml_branch \
            and (raxml_parent.parent.contained_seqids \
              == satchmo_node.contained_seqids \
              or raxml_parent.parent.contained_seqids \
              == all_seqids - satchmo_node.contained_seqids):
          raxml_node = raxml_parent.parent
          match_satchmo_branch_with_raxml_branch(satchmo_node, raxml_node)
          continue
        found_raxml_node = False
        for child in raxml_parent.children:
          if child not in satchmo_branch_of_raxml_branch and \
              (child.contained_seqids == satchmo_node.contained_seqids \
              or child.contained_seqids \
                == all_seqids - satchmo_node.contained_seqids):
            raxml_node = child
            match_satchmo_branch_with_raxml_branch(satchmo_node, raxml_node)
            found_raxml_node = True
            break
        if found_raxml_node:
          continue
      for satchmo_child in satchmo_node.children:
        if satchmo_child in raxml_branch_of_satchmo_branch:
          raxml_child = raxml_branch_of_satchmo_branch[satchmo_child]
          if raxml_child.parent \
              and raxml_child.parent not in satchmo_branch_of_raxml_branch \
              and (raxml_child.parent.contained_seqids \
                == satchmo_node.contained_seqids \
                or raxml_child.parent.contained_seqids \
                == all_seqids - satchmo_node.contained_seqids):
            raxml_node = raxml_child.parent
            match_satchmo_branch_with_raxml_branch(satchmo_node, raxml_node)
            break
          found_raxml_node = False
          for child in raxml_child.children:
            if child not in satchmo_branch_of_raxml_branch and \
                (child.contained_seqids == satchmo_node.contained_seqids \
                or child.contained_seqids \
                  == all_seqids - satchmo_node.contained_seqids):
              raxml_node = child
              match_satchmo_branch_with_raxml_branch(satchmo_node, raxml_node)
              found_raxml_node = True
              break
          if found_raxml_node:
            break
  # First, the SATCHMO root subdivides a branch of the RAxML tree.  At this
  # point one of the two children of the SATCHMO root has a branch length
  # assigned, corresponding to the length of the subdivided branch, and the
  # other has no branch length assigned.  Instead, assign half the length of
  # the subdivided branch to each of these two children.
  satchmo_node = None
  other_node = None
  for child in satchmo_root.children:
    if child in raxml_branch_of_satchmo_branch:
      other_node = child
    else:
      satchmo_node = child
  satchmo_node.branch_length = other_node.branch_length / 2
  other_node.branch_length = other_node.branch_length / 2
  raxml_node = raxml_branch_of_satchmo_branch[other_node].parent
  # At this point, the branches in the SATCHMO tree without lengths are the
  # ones on the path between the root of the SATCHMO tree and the trifurcation
  # of the RAxML tree
  # We will go along this path assigning branch lengths.

  found_unassigned_child = True
  while found_unassigned_child and raxml_node:
    found_unassigned_child = False
    for satchmo_child in satchmo_node.children:
      if satchmo_child not in raxml_branch_of_satchmo_branch:
        if satchmo_child.contained_seqids == all_seqids - raxml_node.contained_seqids:
          found_unassigned_child = True
          match_satchmo_branch_with_raxml_branch(satchmo_child, raxml_node)
          break
    if found_unassigned_child:
      raxml_node = raxml_node.parent
      satchmo_node = satchmo_child

def main():
  if len(sys.argv) < 3:
    print "Usage: %s <rooted_Newick_tree_file> " \
          % sys.argv[0] + " <unrooted_Newick_tree_file_with_lengths>"
    sys.exit(0)

  trivial_translation = string.maketrans('', '')
  f = open(sys.argv[1])
  satchmo_tree_str = f.read().translate(trivial_translation, string.whitespace)
  f.close()
  f = open(sys.argv[2])
  raxml_tree_str = f.read().translate(trivial_translation, string.whitespace)
  f.close()
  satchmo_root = node()
  satchmo_root.readFromTreeString(satchmo_tree_str, 0)
  raxml_root = node()
  raxml_root.readFromTreeString(raxml_tree_str, 0)

  incorporate_branch_lengths_from_unrooted_tree(satchmo_root, raxml_root)

  satchmo_root.printSelf()
  sys.stdout.write("\n")

if __name__ == '__main__':
  main()
