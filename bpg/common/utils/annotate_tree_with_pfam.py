#!/usr/bin/python
# File: annotate_phogs_with_pfam.py
# Author: Grant Shoffner
# Desc: Finds the PHOGs within a given tree and uses hmmscan to detect pfam
#	domains in the sequences in each PHOG. Returns the tree with PHOG-
#	leaves annotated with the domain architecture.


import sys, os
import unittest
import subprocess as sp
from optparse import OptionParser
from Bio import Phylo
from pfacts003.phylofacts.models import Family, Tree, TreeNode, SequenceHeader


# Setting up command line option parsing and error warnings.
Parser = OptionParser()
Parser.add_option('-b', '--book', dest = 'bookname', help = "bpg accession of \
the PF family to annotate, ie bpg0164199")

#Parser.add_option('-m', '--method', dest = 'treemethod', help = "tree type to \
#annotate, either nj or ml")

Parser.add_option('-t', '--tree', dest = 'treefile', help = "tree file to \
annotate, in newick format")

Parser.add_option('-o', dest = 'outtreename', help = "name of file to write \
final tree too", default = False)

Parser.add_option('-c', '--collapse', dest = 'collapsetree', help = "collapse \
clades with common domain structure, write tree to FILE", default = False)

Parser.add_option('--asciitree', dest = 'asciitree', help = "plot tree to \
ascii file", default = False)

Parser.add_option('--names', action = 'store_false', dest = 'newnames', \
default = True, help = "do NOT convert SEQHRD names to UniProt accessions")

Parser.add_option('-w', '--writetbl', action = 'store_true', dest = 'writetable', \
default = False, help = "write Pfams domains for each Phog to a table")

Parser.add_option('-q', '--quiet', action = 'store_false', dest = 'verbose', \
default = True, help = "suppress progress reports")

Parser.add_option('--test', action = 'store_true', dest = 'test', \
default = False, help = "run in testing mode")

(CmdLineOps, Args) = Parser.parse_args()
Test = CmdLineOps.test
Verbose = CmdLineOps.verbose

class TestPhogParsing(unittest.TestCase):
	def setUp(self):
		self.TestFam = Family.objects.get(pk = '0164199')
		self.TestTree = Tree.objects.get(family = self.TestFam, \
							method = 'ml')
		self.TestPhog = TreeNode.objects.get(id = '162118568')
		self.TestUniProtSeqs = [leaf for leaf \
			in self.TestPhog.get_included_leaves() \
                        if leaf.sequence_header.uniprot\
                        and leaf.sequence_header.uniprot.fasta] 
		self.TestHmmscanFile = '/clusterfs/ohana/bpg/\
domain_architecture/testing/PHOG011070253.hmmscan'
		self.TestSeqDKeys = ['SEQHDR20987', 'SEQHDR20038', \
					'SEQHDR20891']
		self.TestSeqDAcc = set(['PF03171.13'])
		self.TestPhyloTree = Phylo.read('/clusterfs/ohana/bpg/\
domain_architecture/testing/bpg0100031.nj', 'newick')
		self.TestTargetLeaves = set(['SEQHDR20323', 'SEQHDR20132', 'SEQHDR20227'])

	def test_correct_number_phogs_retrieved(self):
		self.assertEqual(len(getphogsfromtree(self.TestTree)), 242)

	def test_getphogsfromtree_does_not_return_dict(self):
		self.assertRaises(TypeError, \
				getphogsfromtree(self.TestTree), {})

	def test_writesequencefile_returns_seqs(self):
		ReturnedSeqs = writesequencefile(self.TestPhog, \
						ReturnSeqs = True)
		self.assertEqual(ReturnedSeqs, self.TestUniProtSeqs)

	def test_writesequencefile_on_example_phog(self):
		try: StdOutputFile = open('/clusterfs/ohana/bpg/\
domain_architecture/testing/162118568.fa')
		except IOError:
			print '\n\tExample sequence file not found for PHOG %s' \
				% self.TestPhog.__str__()
			return
		PhogFileName = self.TestPhog.__str__() + '.fa'
		if PhogFileName in os.listdir(os.getcwd()):
			print '\n\tPHOG sequence file already exists in \
working directory'
			return
		writesequencefile(self.TestPhog)
		PhogSeqFile = open(PhogFileName)
		PhogSeqFileLines = PhogSeqFile.readlines()
		PhogSeqFile.close()
		StdOutputFileLines = StdOutputFile.readlines()
		StdOutputFile.close()
		os.remove(PhogFileName)
		self.assertEqual(PhogSeqFileLines, StdOutputFileLines)

	def test_convertseqhdr_works_on_example(self):
		self.assertEqual(convertseqhdr('SEQHDR17540') ,'O00469')

	def test_parsehmmscanoutput_returns_correct_keys(self):
		SeqDomains = parsehmmscanoutput(self.TestHmmscanFile)
		self.assertEqual(SeqDomains.keys(), self.TestSeqDKeys)

	def test_parsehmmscanoutput_finds_correct_acc(self):
		SeqDomains = parsehmmscanoutput(self.TestHmmscanFile)
		JustDomainAccs = set()
		for Seq in SeqDomains.keys():
			for Domain in SeqDomains[Seq]:
				JustDomainAccs.add(Domain[0])
		self.assertEqual(JustDomainAccs, self.TestSeqDAcc)

	def test_parsehmmscanoutput_correct_length(self):
		SeqDomains = parsehmmscanoutput(self.TestHmmscanFile)
		for Seq in SeqDomains.keys():
			for Domain in SeqDomains[Seq]:
				self.assertEqual(len(Domain), 4)

	def test_getparent_works_on_leaf(self):
		leaf = self.TestPhyloTree.get_terminals()[0]
		parent = getparent(leaf, self.TestPhyloTree)
		ParentLeaves = set([node.name for node \
					in parent.get_terminals()])
		self.assertEqual(ParentLeaves, self.TestTargetLeaves)

	def test_getparent_works_at_root(self):
		path = self.TestPhyloTree.get_path(self.TestTargetLeaves.pop())
		node = path[0]
		self.assertEqual(getparent(node, self.TestPhyloTree), \
						self.TestPhyloTree.root)

TestingSuite = unittest.TestLoader().loadTestsFromTestCase(TestPhogParsing)

def getphogsfromtree(TreeObject):
	'''Retrieve the PHOG-Ss associated with a given tree.'''
	# BROKEN: Needs to be updated with Ruchira's new code for getting PHOGs
	# --GS 22.3.11
	InitialNodes = TreeNode.objects.filter(tree = TreeObject.id, \
				superorthologous_node__isnull = False)
	NodeIds = set([treenode.superorthologous_node.id for treenode \
				in InitialNodes])
	FilteredNodes = set([TreeNode.objects.get(id = ID) for ID in NodeIds])
	return FilteredNodes

def writesequencefile(treenode, ReturnSeqs = False):
	'''Takes a TreeNode object and writes the full sequences in the\
tree to a file.'''
	TreeName = treenode.__str__()
	UniProtSeqs = [leaf for leaf in treenode.get_included_leaves() \
			if leaf.sequence_header.uniprot\
			and leaf.sequence_header.uniprot.fasta]
	# Current naming convention for the sequence files is 
	# "<Phog_Accession>.fa". If this is changed, make sure to change 
	# RunHmmscan also.
	SeqFile = open(TreeName + '.fa', 'a')
	for Seq in UniProtSeqs:
		SeqFile.write('>SEQHDR' + repr(Seq.sequence_header.id) + '\n')
		SeqFile.write(Seq.sequence_header.uniprot.chars + '\n')
	SeqFile.close()
	if ReturnSeqs: return UniProtSeqs
	else: return

def runhmmscan(treenode):
	'''Starts an hmmscan subprocess agianst Pfam-A on the sequence file \
for the given TreeNode. Returns the name of the file the results are written to.'''
	# Note: this function requires (1) that writesequencefile(treenode) be run
	# first, and (2) that the name of the sequence file for treenode follows
	# the convention "<treenode.__str__()>.fa".
	SeqFileName = treenode.__str__() + '.fa'
	TblFileName = treenode.__str__() + '.hmmscan'
	HmmscanProcessName  = 'hmmscan --domtblout %s -o /dev/null -E 0.001 \
--domE 0.001 /clusterfs/ohana/external/pfam/current/Pfam-A.hmm %s' % \
		(TblFileName, SeqFileName)
	sp.Popen(HmmscanProcessName, shell = True).wait()
	return TblFileName

def parsehmmscanoutput(OutputFileName):
	'''Parses the table output of hmmscan, returns a dictionary with\
sequence ids as keys to lists of the domains in the sequence.'''
	# Dictionary keys gives tuples (Domain_Accession, Eval, Start, End)
        OutputFile = open(OutputFileName)
        SeqDomains = {}
        for Line in OutputFile.readlines()[3:]:
                Line = Line.split()
                SeqName = Line[3] # Sequence name in SEQHDR format
                DomainAcc = Line[1] # Pfam accession
		DomainIE = float(Line[12]) # Independant E-value for this domain
		if DomainIE > 0.001: continue
                DomainStart = int(Line[17]) # Start index of hit
                DomainEnd = int(Line[18]) # End index of hit
                Domain = (DomainAcc, DomainIE, DomainStart, DomainEnd)
                if SeqName in SeqDomains.keys():
                        SeqDomains[SeqName].append(Domain)
                else: SeqDomains[SeqName] = [Domain]
	OutputFile.close()
        for Seq in SeqDomains.keys():
		ToBeRemoved = []
		for Domain in SeqDomains[Seq]:
			if Domain in ToBeRemoved: continue
			DStart = Domain[2]
			DEnd = Domain[3]
			DEval = Domain[1]
			for TDomain in SeqDomains[Seq]:
				if TDomain in ToBeRemoved: continue
				if Domain == TDomain: continue
				if (TDomain[2] <= DStart <= TDomain[3]\
					or TDomain[2] <= DEnd <= TDomain[3])\
					and TDomain[1] > DEval:
						ToBeRemoved.append(TDomain)
		SeqDomains[Seq] = sorted([Domain for Domain in SeqDomains[Seq] \
					if Domain not in ToBeRemoved], \
					key = lambda DomList: DomList[2])
        return SeqDomains

def writephogdomains(PhogName, SeqDomains):
	TableFile = open(PhogName + '_domain.tbl', 'a')
	for SeqName in SeqDomains.keys():
		TableFile.write(SeqName + ' ')
		for Domain in SeqDomains[SeqName]:
			TableFile.write(Domain[0] + ' ')
		TableFile.write('\n')
	TableFile.close()
	return

def convertseqhdr(SeqHdrName):
	'''Converts a leaf label from the SEQHRD###### format to the UniProt\
accession, if possible.'''
	SeqHdrId = SeqHdrName[6:]
	Leaf = SequenceHeader.objects.get(id = SeqHdrId)
	if Leaf.uniprot: ReturnName = Leaf.uniprot.accession
	else: ReturnName = SeqHdrName
	return ReturnName

def rewritetreeleaveswithdomains(PhyloTree, SeqDomains):
	'''Takes a parsed tree, a dictionary of SEQHRDs and domains,\
and returns a tree with the leaves relabeled with the domains.'''
	WorkingTree = PhyloTree
	TerminalNodes = WorkingTree.get_terminals()
	for i in xrange(len(TerminalNodes)):
		NodeName = TerminalNodes[i].name
		if NodeName in SeqDomains.keys():
			NewName = NodeName + ' - ' + ' '.join([Domain[0] for \
						Domain in SeqDomains[NodeName]])
			WorkingTree.get_terminals()[i].name = NewName
	return WorkingTree

def getparent(Node, PhyloTree):
	'''Finds the parent node of a node in a tree.'''
	path_to_root = PhyloTree.get_path(Node)
	if len(path_to_root) == 1: Parent = PhyloTree.root
	else: Parent = path_to_root[-2]
	return Parent

def findmaximalsubtrees(PhyloTree, SeqDomains):
	'''Finds the maximal subtrees of a tree for which all terminal nodes\
share a common domain structure.'''
	# Define a dictionary with just the domain names, not the evalues or 
	# start/end indices.
	JustSeqDomains = {}
	for SeqName in SeqDomains.keys():
		JustSeqDomains[SeqName] = tuple([Domain[0] for Domain \
						in SeqDomains[SeqName]])
	LeaveswithUniProt = set([leaf for leaf in PhyloTree.get_terminals() \
				if leaf.name in JustSeqDomains.keys()])
	CoveredLeaves = set()
	MaximalNodes = set()
	for leaf in LeaveswithUniProt:
		if leaf in CoveredLeaves: continue
		maximalnode = leaf
		currentnode = leaf
		while True:
			parent = getparent(currentnode, PhyloTree)
			children = set(parent.get_terminals())\
.intersection(LeaveswithUniProt)
			Domains = set([JustSeqDomains[node.name] for node in children])
			if len(Domains) != 1:
				MaximalNodes.add(maximalnode)
				for node in maximalnode.get_terminals():
					CoveredLeaves.add(node)
				break
			if parent == PhyloTree.root:
				MaximalNodes.add(parent)
				for node in parent.get_terminals():
					CoveredLeaves.add(node)
				break
			maximalnode = parent
			currentnode = parent
	return MaximalNodes

def collapsemaximalsubtrees(PhyloTree, MaximalNodes, SeqDomains):
	WorkingTree = PhyloTree
	JustSeqDomains = {}
	for SeqName in SeqDomains.keys():
		JustSeqDomains[SeqName] = tuple([Domain[0] for Domain \
						in SeqDomains[SeqName]])
	for maximalnode in MaximalNodes:
		ConsensusDomains = set([JustSeqDomains[leaf.name] for leaf \
					in maximalnode.get_terminals() \
					if leaf.name in JustSeqDomains.keys()])
		NodeLabel = ' '.join(ConsensusDomains.pop())
		if maximalnode in PhyloTree.get_terminals():
			maximalnode.name += ' - ' + NodeLabel
			continue
		maximalnode.collapse_all()
		for leaf in maximalnode.get_terminals():
			maximalnode.collapse(leaf)
		maximalnode.name = NodeLabel
	return WorkingTree



if __name__ == '__main__':
	# Checking for proper user input.
	if CmdLineOps.bookname is None and not Test:
		print 'Missing mandatory PF book accession argument.'
		sys.exit(2)
	#if CmdLineOps.treemethod not in ('ml','nj') and not Test:
        #	print 'Missing mandatory tree type (nj or ml) argument.'
        #	sys.exit(2)
	if CmdLineOps.treefile is None and not Test:
		print 'Missing mandatory tree file name argument.'
		sys.exit(2)

	# Run unittests if given command line option.
	if Test:
		if Verbose: print 'Starting in testing mode.'
		unittest.TextTestRunner(verbosity=2).run(TestingSuite)
		sys.exit(0)

	# Retrieve the desired PhyloFacts family
	try: Fam = Family.objects.get(pk = CmdLineOps.bookname[3:])
	except: 
		print 'Family not found. Is your bpg accession formated correctly?'
		sys.exit(1)

        # Open the newick tree file to annotate.
	try: TreetoAnnotate = Phylo.read(CmdLineOps.treefile, 'newick')
	except IOError:
		print 'Tree file not found.'
		sys.exit(1)

	# Write the sequences in the family to a file and run hmmscan
	if Verbose:
		sys.stdout.write('Writing sequences for Family %s to file ... '\
				% Fam.get_accession())
		sys.stdout.flush()

	writesequencefile(Fam.canonical_root_node())

	if Verbose: 
		print 'done'
		sys.stdout.write('Running hmmscan against Pfam-A ... ')
		sys.stdout.flush()

	hmmtblfilename = runhmmscan(Fam.canonical_root_node())

	if Verbose:
		print 'done'
		sys.stdout.write('Parsing hmmscan output ... ')
		sys.stdout.flush()

	SequenceDomains = parsehmmscanoutput(hmmtblfilename)

	if Verbose: print 'done'	

	# Write a Sequence <--> Domain association table if given the option
	if CmdLineOps.writetable:
		if Verbose:
			sys.stdout.write('Writing Pfam domains to table file ... ')
			sys.stdout.flush()

		writephogdomains(Fam.canonical_root_node().__str__(), 
							SequenceDomains)

		if Verbose: print 'done'

	# Rewrite the newick tree with new leaf labels and print an ascii tree
	if Verbose:
		sys.stdout.write('Annotating tree with Pfam domains ... ')
		sys.stdout.flush()

	if CmdLineOps.outtreename:
		AnnotatedTree = rewritetreeleaveswithdomains(TreetoAnnotate, 
							SequenceDomains)
		Phylo.write(AnnotatedTree, CmdLineOps.outtreename, 'newick',
						branchlengths_only=True)

	if CmdLineOps.asciitree:
		OutTreeFile = open(CmdLineOps.outtreename, 'w')
		Phylo.draw_ascii(AnnotatedTree, OutTreeFile)
		OutTreeFile.close()

	if Verbose: 
		if CmdLineOps.outtreename or CmdLineOps.asciitree: print 'done'
		else: print 'skipped'

	if CmdLineOps.collapsetree:
		# Find the maximal subtrees with consensus domain structure
		if Verbose:
			sys.stdout.write('Collapsing subtrees with common \
domain structure ... ')
			sys.stdout.flush()
		MaximalNodes = findmaximalsubtrees(TreetoAnnotate, SequenceDomains)
		CollapsedTree = collapsemaximalsubtrees(TreetoAnnotate, \
							MaximalNodes, \
							SequenceDomains)
		Phylo.write(CollapsedTree, CmdLineOps.collapsetree, 'newick',
						branchlengths_only=True)
		if Verbose: print 'done'


        #----------------------------
        #TerminalNames = [node.name for node in TreetoAnnotate.get_terminals()]
        #print TerminalNames == SequenceDomains.keys()
        #print '-------------'
        #print TerminalNames
        #print '-------------'
        #print SequenceDomains.keys()
        #sys.exit(1)
        #----------------------------

	#----------------------------
	# Some testing code!
	#for max_node in MaximalNodes:
	#	print '---------------------------------------------------'
	#	for leaf in max_node.get_terminals():
	#		print leaf.name
	#----------------------------
