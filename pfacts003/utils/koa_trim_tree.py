#from ete2 import Tree
#from ete2 import PhyloTree

from __future__ import division
import csv
import itertools
import functools
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment

class MSATree:
    """An MSA with a Phylo Tree"""
    def __init__(self, msa, tree):
        self.msa = msa
        self.tree = tree

        # We pull what is before the white space to match with the tree
        def get_key(align):
            return align.name.split(' ')[0]

        self.msa_by_name = dict((get_key(alignment), alignment)
                                for alignment in self.msa)

        # stolen from http://biopython.org/wiki/Phylo_cookbook
        self.parents = {}
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                self.parents[child] = clade

        # Flag for if we add the root to the output trees
        self.added_root = False

    def _node_distance(self, first, second):
        """Takes in two clade nodes

        Returns the pairwise identity percentage

        Nodes are remembered by name, and the results cached
        """

        name_1 = first.name.split(' ')[0]
        name_2 = second.name.split(' ')[0]

        seq1 = self.msa_by_name[name_1]
        seq2 = self.msa_by_name[name_2]

        distance = self._seq_distance(seq1, seq2)

        return distance



    class memoized(object):
       def __init__(self, func):
          self.func = func
          self.cache = {}
       def __call__(self, *args):
          try:
             return self.cache[args]
          except KeyError:
             value = self.func(*args)
             self.cache[args] = value
             return value
          except TypeError:
             return self.func(*args)
       def __repr__(self):
          return self.func.__doc__
       def __get__(self, obj, objtype):
          return functools.partial(self.__call__, obj)
    
    @memoized
    def _seq_distance(self, seq1, seq2):
        num_positions = 0  # the number of positions with at least one amino
        num_matches = 0    # the number of matches seen

        for pair in zip(seq1, seq2):
            if pair[0].islower() or pair[1].islower() \
                    or pair[0] == "." or pair[1] == "." \
                    or ( pair[0] == "-" and pair[1] == "-" ):
                continue
    
            else:
                num_positions += 1
    
                if (pair[0] == pair[1]):
                    num_matches += 1
    
        return num_matches / float(num_positions) * 100


    def subtree_distances(self, root):
        """Takes a node

        returns a list of the percent identities between all leaves

        This isn't needed for the program to run, but can help with debugging
        The other method is faster, as it only computes distances until one
        exeeds the threshold
        """

        nodes = root.get_terminals()
        nodes.reverse()
        node_pairs = itertools.ifilter(
            lambda (a1, a2): a1.name < a2.name,
            itertools.product(nodes, nodes))

        distances = [self._node_distance(pair[0], pair[1])
            for pair in node_pairs]

        return distances

    def _subtree_above_threshold(self, root, threshold=0):
        """Takes in a clade node and a threshold

        Returns True if all children have a pairwise sequence
        alignment below the threshold"""

        # BFS will grab the most distant leaves last
        # reverse should guarantee that we get the smallest subtrees
        nodes = root.get_terminals()
        nodes.reverse()

        # always making n1 < n2 so that we can memoize
        # this all could be done in parallel
        node_pairs = itertools.ifilter(
            lambda (a1, a2): a1.name < a2.name,
            itertools.product(nodes, nodes))

        for pair in node_pairs:
            distance = self._node_distance(pair[0], pair[1])
            if distance and distance < threshold:
                return False

        return True

    def _subtree_below_maximum_leaves(self, root, threshold):
        """ takes a clade node and a threshold
        
        Returns true if there are fewer or equal leaves than threshold
        Returns false otherwise """

        nodes = root.get_terminals()
        return len(nodes) <= threshold


    def _meets_criteria(self, root, criteria):
        if not self._subtree_below_maximum_leaves(root, criteria):
            return False

        return True


    def maximal_subtree_msa(self, start, criteria):
        for node in self.tree.get_terminals():
            if node.name.split(' ')[0] == start:
                starting_node = node
                break

        if not starting_node:
            raise Exception("Cannot find seed node in tree")

        node = starting_node

        while self.parents.get(node) and self._meets_criteria(self.parents.get(node), criteria):
            node = self.parents.get(node)

        new_tree = Phylo.BaseTree.Tree.from_clade(node)

        def get_alignment_from(tree):
            """Fetches the MSA for a Phylo Tree"""
            msa = []
            for node in tree.get_terminals():
                alignment = self.msa_by_name[node.name.split(' ')[0]]
                if msa:
                    msa.append(alignment)
                else:
                    msa = MultipleSeqAlignment([alignment])

            return msa

        return [ new_tree, get_alignment_from(new_tree) ]


def main():
#    unrooted_tree = open('tree.newick', 'r').read()
#    t = PhyloTree("tree.newick", alignment="chopped.fa", alg_format="fasta")
#	test_seqs = AlignIO.read("chopped.fa", "fasta")
#	print fast_tree(test_seqs)
    msa_tree = MSATree(
        AlignIO.read("hmmer_top500_mafft.fa", "fasta"),
        Phylo.read("hmmer.tree", "newick"))

    [ subtree, subtree_msa ] = msa_tree.maximal_subtree_msa("tr|B2CP06|B2CP06_MOUSE", 100)
    print subtree.format("newick")
        
if __name__ == "__main__":
	main()