#!/usr/bin/env python2.6

help_str = '''
Script takes as input a tree, an alignment used to construct the tree (headers should be present)
in the tree, the % SH test threshold, and the % min pwid

e.g. ../bpg0240116.tree ../bpg0240116.a2m 0.7 1.0

and prints out subtrees after doing the following
* Parse the tree and the alignment.
* Associate all the tree nodes with the corresponding records in the alignment.
* For every tree node, store the SH test scores and the minimum pairwise identity of the records contained within that node.
* Traverse down from the root and stop at internal nodes resulting in the largest subtrees with minimum pairwise identity >= %pwid.
* For each of these nodes, traverse up to the root and stop at the last node where the sh test score < %sh
* Finally, for each child of the resulting nodes, print
** A node id (which could be tree_node1...n)
** All the sequence headers contained within that node.Parse the tree and the alignment.

'''
import os
import sys
import tempfile
from ete2 import PhyloTree
from StringIO import StringIO
sys.path.insert(0, '/home/awarrier/ohana_repository')
from bpg.common.BPGPWID import pairwise_identity_belvu as pwid
from itertools import combinations

def integrate_pwids_into_tree(tree, alignment):
    '''Takes a tree and an alignment and returns a new tree with values of pwid added to each
    node in the tree as node.pwid.'''
    pt = PhyloTree(tree,alignment=alignment,alg_format="fasta")
    for ind, node in enumerate(pt.traverse()):
        node.node_kerf_name = 'node%s' % str(ind).zfill(3)
        # For later kerf and sh functions
        node.kerf_pass = False
        node.sh_pass = False
        if not node.is_leaf():
            node.min_pwid = get_min_pwid_of_leaves(node.get_leaves())
        else:
            node.min_pwid = 1.0
    return pt
        
def get_min_pwid_of_leaves(leaves):
    '''Takes leaves and returns their minimum_pairwise_identity_belvu
    '''
    leaf_sequences = [leaf.sequence for leaf in leaves]
    pairwise_identity_list = [pwid(seq0, seq1) for seq0, seq1
                              in combinations(leaf_sequences, 2)]
    return min(pairwise_identity_list)

def clean_alignment(alignment_file):
    '''Sometimes the alignment needs to be cleaned so that the headers match the leaf
    labels. All code for that goes in here.
    Input arguments: path to alignment file
    output: path to temporary cleaned alignment file.
    '''
    f = open(alignment_file).read().splitlines()
    output_list = []
    for line in f:
        if line.startswith('>'):
            output_list.append(line.split()[0])
        else:
            output_list.append(line)
    outfile = tempfile.NamedTemporaryFile(mode='w+b')
    outfile.write('\n'.join(output_list))
    outfile.seek(0)
    return outfile

def pwid_match(node):
    if node.min_pwid >= 70:
        return True
    else:
        return False

def kerf(node, pwid_thr):
    '''Takes a tree and does the first part of kerf, i.e. stores True where the internal nodes
    are kerf supported.
    with minimum_pairwise_identity greater than or equal to % min pwid argument.
    '''
    if not node.is_leaf():
        if node.min_pwid >= pwid_thr:
            return node
    for child in node.children:
        kerf_node = kerf(child, pwid_thr)
        if kerf_node:
            kerf_node.kerf_pass = True
            continue
def sh(tree, sh_thr):
    '''Takes a tree, sets the sf_pass flag to True only if its own support is less than threshold,
    if either child node has an sh support value greater than threshold and if the child with an
    sh support value greater than threshold contains a node with kerf_pass set to True.
    '''
    for node in tree.traverse():
        if node.support < sh_thr:
            for child in node.children:
                if child.support > sh_thr:
                    kerf_pass_list = [desc.kerf_pass for desc in child.get_descendents()
                                      if desc.ker_pass]
                    if kerf_pass_list:
                        node.sh_pass = True

def traverse_to_sh_nodes(nodes):
    '''Takes a list of node objects, generated by kerf(tree) and traverses up the tree
    until the sh test value is less than the % SH test threshold argument.
    Then return the child nodes for these.
    '''
    pass

def print_sh_nodes(tree):
    '''Prints the result out.
    '''
    for node in tree.traverse():
        if node.sh_pass:
            for child in node.children:
                print ';;', child.name
                for leaf in child.get_leaves():
                    print leaf.name
                

def main(tree, alignment, sh_thr_pct, min_pwid_thr_pct):
    '''
    '''
    new_alignment = clean_alignment(alignment)
    pwid_tree = integrate_pwids_into_tree(tree, new_alignment.name)
    new_alignment.close()
    kerf(pwid_tree, min_pwid_thr_pct)
    sh(pwid_tree, sh_thr_pct)
    print_sh_nodes(pwid_tree)

if __name__ == "__main__":
    try:
        tree, alignment, sh_thr_pct, kerf_thr = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    except:
        print """All four arguments are needed.
        Argument1 = path to newick tree
        Argument2 = path to alignment created from this tree
        Argument3 = sh threshold
        Argument4 = kerf threshold value.\n%s""" % help_str
        sys.exit()
        
    main(tree, alignment, float(sh_thr_pct), float(kerf_thr))
