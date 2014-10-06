#!/usr/bin/python
# File: redundancyfilter.py
# Desc: Reduce the redundancy of sequences in a group of PhyloFacts families 
#	sharing a common Pfam domain.


import sys, os, csv, re, logging, threading, Queue
from optparse import OptionParser
import Bio.SeqIO
import pfacts003.phylofacts.models as pfacts


class ParseMSA(threading.Thread):
	"""Parses the identifier for each sequence in an MSA."""
	msaQueue = Queue.Queue()
	lock = threading.Lock()
	accessions = {}

	def __init__(self):
		self.threadid = os.urandom(8)
		threading.Thread.__init__(self)

	def run(self):
		while True:
			try:
				MsaID, MsaFileName = ParseMSA.msaQueue.get(block = False)
			except Queue.Empty:
				break
			MsaRecords = Bio.SeqIO.parse(open(MsaFileName), 'fasta')
			#SequenceIdentifiers = [record.id.split('|')[1] \
			#			for record in MsaRecords]
			SequenceIdentifiers = []
			for record in MsaRecords:
				try: accession = record.id.split('|')[1]
				except IndexError:
					accession = record.id.split('_')[0]
				SequenceIdentifiers.append(accession)


			if not SequenceIdentifiers:
				logging.warning('No sequences parsed for %s' \
									% MsaID)
			ParseMSA.lock.acquire()
			ParseMSA.accessions[MsaID] = SequenceIdentifiers
			ParseMSA.lock.release()
		return


def parse_input_file(csvfilename):
	if not os.path.isfile(csvfilename):
		logging.error('Input file %s not found.' % csvfilename)
		sys.exit(1)

	PfamAccession = re.compile('PF\d\d\d\d\d')
	DomainFamilies = {}
	csvfile = csv.reader(open(csvfilename))
	linecounter = 0

	for row in csvfile:
		linecounter += 1
		if not row:
			logging.info('Empty row near line %g in input %s'
					% (linecounter, csvfilename))
			continue
		accession = row[0]
		if len(accession) != 7 or not PfamAccession.match(accession):
			logging.warning(
'Bad Pfam accession near line %g, skipping' % linecounter)
			logging.info(
'Bad Pfam accession: %s, line %g, in input %s' % \
	(accession, linecounter, csvfilename))
			continue
		if len(row) == 1:
			if accession not in DomainFamilies:
				DomainFamilies[accession] = []
		elif len(row) == 2:
			MsaFile = row[1].strip()
			if not os.path.isfile(MsaFile):
				logging.error('MSA file %s not found' 
						% MsaFile)
				sys.exit(1)
			if accession not in DomainFamilies.keys():
				DomainFamilies[accession] = [MsaFile]
			else:
				DomainFamilies[accession].append(MsaFile)

		elif len(row) > 2:
			logging.warning('Incorrect number of columns near \
line %g in input file, skipping this entry.' % linecounter)
			continue
	return DomainFamilies

def get_families_with_domain(pfam_domain_accession, return_acc = False):
	"""Retrieves family objects from PF3.0 containing the given domain"""
	try: pfam_domain = pfacts.Pfam.objects.get(
		accession = pfam_domain_accession, overall_pfam_version=24)
	except pfacts.Pfam.DoesNotExist:
		logging.error(
'No Pfam object found in database for accession %s' % pfam_domain_accession
				)
		raise
	try: pfam_hmm = pfacts.HMM.objects.get(pfam = pfam_domain)
	except pfacts.HMM.DoesNotExist:
		logging.error(
'No HMM object found for Pfam accession %s' % pfam_domain_accession
				)
		raise
	sequences_with_domain = pfacts.SequenceHMM.objects.filter(
								hmm = pfam_hmm)
	tree_node_consensus = pfacts.TreeNodeConsensus.objects.filter(
		sequence__in = [hit.sequence for hit in sequences_with_domain])
	families_with_domain = [consensus.tree_node.tree.family for consensus \
			in tree_node_consensus \
			if consensus.tree_node.tree.family.status != 'bad' and consensus.tree_node.tree.family.family_type_id == 'C']
	if not families_with_domain:
		logging.info('Found no families for domain %s' \
						% pfam_domain_accession)
	if return_acc:
		return [fam.get_accession() for fam in families_with_domain]
	else:
		return families_with_domain

def get_family_msa(accession):
	msafile = '/clusterfs/ohana/bpg/pfacts/' + '/'.join([accession[:4], 
							accession[:7], 
							accession, 
							accession + '.a2m'])
	if not os.path.isfile(msafile):
		logging.error(
'Could not find MSA file for family %s and path %s' % (accession, msafile))
		sys.exit(1)
	return msafile


def get_sequence_accessions(family):
	#Currently not needed, but might be handy...
	sequences = family.canonical_root_node().get_included_leaves()
	sequences = [seq for seq in sequences \
			if seq.sequence_header.uniprot \
			and seq.sequence_header.uniprot.accession]
	accessions = set([seq.sequence_header.uniprot.accession \
						for seq in sequences])
	return accessions


def filterfamilies(families, threshold):
	"""Returns the smallest subset of families covering these sequences."""
	#Start parsing of family MSA files
	for family in families:
		ParseMSA.msaQueue.put((family, family))
	Threads = [ParseMSA() for _ in range(3)]
	for thread in Threads:
		thread.start()
	for thread in Threads:
		thread.join()

	#Aggregate all the UniProt accessions
	accessions = []
	fam_acc = {}
	for family in families:
		accessions_for_this_fam = ParseMSA.accessions[family]
		accessions += accessions_for_this_fam
		fam_acc[family] = set(accessions_for_this_fam)
	accessions = set(accessions)

	#Find the covering of this set of accessions
	Covering = []
	CoverageStats = []
	families_remaining = list(families)
	accessions_to_cover = set(accessions)
	while len(accessions_to_cover) > 0:
		maximum_family = ''
		maximum_coverage = 0
		for family in families_remaining:
			coverage = len(accessions_to_cover & fam_acc[family])
			if coverage > maximum_coverage:
				maximum_family = family
				maximum_coverage = coverage
		if maximum_family and maximum_coverage >= threshold:
			Covering.append(maximum_family)
			CoverageStats.append((maximum_family, 
					repr(len(fam_acc[maximum_family])), 
					repr(maximum_coverage)))
			families_remaining.remove(maximum_family)
			accessions_to_cover = accessions_to_cover - \
						fam_acc[maximum_family]
		else:
			break
	return Covering, CoverageStats, fam_acc

def print_coverage_data(outputfile, DomainAccessions, Covering):
	total_domains = len(DomainAccessions.keys())
	outputfile.write('---Coverage Data---\n')
	outputfile.write('Total domains filtered: %d' % total_domains)
	for domain in DomainAccessions:
		accessions = []
		for family in DomainAccessions[domain]:
			accessions += ParseMSA.accessions[family]
		accessions = set(accessions)
		accessions_remaining = set(accessions)
		for family in Covering:
			accessions_remaining = accessions_remaining - \
					DomainAccessions[domain][family]
		
		outputfile.write('%s ---\n' % domain)
		outputfile.write('\tTotal unique accessions: %d\n' % \
							len(accessions))
		outputfile.write('\tAccessions left uncovered: %d\n' % \
						len(accessions_remaining))
	return
		
def main():
	Parser = OptionParser('usage: redundancyfilter.py [options] <domains_file>')
	Parser.add_option('-t', '--threshold',
		dest = 'threshold',
		type = 'int',
		default = 10,
		help = "Max new sequences a family must contribute")
	Parser.add_option('-q', '--quiet',
		action = 'store_false',
		dest = 'verbose', 
		default = True, 
		help = "suppress progress reports")
	Parser.add_option('-o', '--output',
		dest = 'output_file',
		default = False,
		help = "send results to FILE, instead of standard output")
	Parser.add_option('-s', '--stats',
		dest = 'stats',
		action = 'store_true',
		default = False,
		help = "print coverage statistics")
	(CmdLineOps, Args) = Parser.parse_args()
	Verbose = CmdLineOps.verbose
	if len(Args) == 1:
		input_file = Args[0]
	else:
		logging.error('Recieved incorrect arguments')
		sys.exit(1)

	if Verbose:
		sys.stdout.write('Parsing input and retrieving exisiting families ... ')
		sys.stdout.flush()
	DomainFamilies = parse_input_file(input_file)
	for domain in DomainFamilies:
		families = get_families_with_domain(domain, return_acc = True)
		family_msas = map(get_family_msa, families)	
		DomainFamilies[domain] += family_msas
	if Verbose:
		print 'done'
		total_domains = len(DomainFamilies.keys())
		total_families = 0
		for domain in DomainFamilies:
			for fam in DomainFamilies[domain]:
				total_families += 1
		print 'Starting with %d domains and a total of %d families' % \
				(total_domains, total_families)
	if CmdLineOps.output_file:
		outputfile = open(CmdLineOps.output_file, 'a')
	else:
		outputfile = sys.stdout

	DomainAccessions = {}
	for domain in DomainFamilies:
		if Verbose:
			sys.stdout.write('Finding covering of domain %s ... ' % domain)
			sys.stdout.flush()

		Covering, CoverageStats, family_accessions = filterfamilies(
							DomainFamilies[domain],
							CmdLineOps.threshold)
		DomainAccessions[domain] = family_accessions
		if Verbose:
			print 'done'

		if CmdLineOps.stats:
			lines_to_print = ['\t'.join(line) for line in \
								CoverageStats]
		else:
			lines_to_print = Covering

		outputfile.write('Domain: %s\n' % domain)
		outputfile.write('Families to retain:\n')
		for family in lines_to_print:
			outputfile.write(family + '\n')
		outputfile.write('Families to discard:\n')
		for family in DomainFamilies[domain]:
			if family not in Covering:
				outputfile.write(family + '\n')

if __name__ == '__main__':
	main()
