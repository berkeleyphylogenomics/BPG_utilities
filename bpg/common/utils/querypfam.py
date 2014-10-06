#!/usr/bin/python
# File: querypfam.py
# Author: Grant Shoffner
# Desc: Defines threaded classes for retrieving and processing sequence
# data from Pfam.

import os, sys, urllib, urllib2, threading, Queue
import xml.dom.minidom
import logging

def pfamdomains(uniprot_accession):
    GetPfamDomains.uidQueue.put(uniprot_accession)
    thread = GetPfamDomains()
    thread.start()
    thread.join()
    return GetPfamDomains.SeqData[uniprot_accession]


class PfamServerError(Exception):
	'''Exception raised when Pfam sequence search returns 
	a bad HTTP code.'''

	def __init__(self, seqid, result_url, threadid, httpcode):
		self.seqid = seqid
		self.result_url = result_url
		self.threadid = repr(threadid)
		self.code = repr(httpcode)

	def __str__(self):
		ErrorMessage = '''Pfam server failed on sequence
Sequence ID: %s
Result URL: %s
HTTP request code: %s
Thread ID: %s

See http://pfam.sanger.ac.uk/help#tabview=tab10 for an explaination of 
the HTTP error codes.'''
		return ErrorMessage % (self.seqid, self.result_url, self.code, \
					self.threadid)

class GetPfamDomains(threading.Thread):
	uidQueue = Queue.Queue()
	SeqData = {}
	seqdata_lock = threading.Lock()

	def __init__(self):
		self.threadid = os.urandom(8)
		threading.Thread.__init__(self)

	def run(self):
		while True:
			try: uid = GetPfamDomains.uidQueue.get(block = False)
			except Queue.Empty:
				break
			url = 'http://pfam.sanger.ac.uk/protein/%s?output=xml' \
				 % uid
			request = urllib2.urlopen(url)
			GetPfamDomains.seqdata_lock.acquire()
			GetPfamDomains.SeqData[uid] = \
				PfamSeqSearch.parsepfammatchxml(request)
			GetPfamDomains.seqdata_lock.release()

class GetUniProtSeq(threading.Thread):
	uidQueue = Queue.Queue()
	Sequences = []
	seqdata_lock = threading.Lock()

	def __init__(self):
		self.threadid = os.urandom(8)
		threading.Thread.__init__(self)

	def run(self):
		while True:
			try: uid = GetUniProtSeq.uidQueue.get(block = False)
			except Queue.Empty:
				break
			url = 'http://www.uniprot.org/uniprot/%s.fasta' \
				% uid
			request = urllib2.urlopen(url)
			GetUniProtSeq.seqdata_lock.acquire()
			GetUniProtSeq.Sequences.append(request)
			GetUniProtSeq.seqdata_lock.release()

TestProtein = ('A3QK15', 
'''MSKDTEAKSEEIMESKVLWYPDSKRNTQTDRFRTLVNREFGLNLANYNDLYQWSVDSYPE
FWAQVWKFCGITCSKMYEEVVDVSKRISDVPEWFKGSRLNYAENLLKHKDQDKVALYAAS
EAKEEIVKVTFGELRRDVALFAAAMRKMGIKIGDRVVGYLPNGVHAVEAMLAAASIGAIW
SSTSPDFGVNGVLDRISQIQPKLIFSVAAVVYNGKQHDHMEKLQNVVKGLPDLKKVVVIP
YVRSRQETDLSKIPNSVFLEDFLATGKEGDQDPQLEFEQLPFSHPLFIMYSSGTTGAPKC
MVHSAGGTLIQHLKEHILHGNMTFNDVIIYYTTTGWMMWNWLISSLAVGASVVLYDGSPL
VPSANVLWDLVDRLGITIFGTGAKWLAVLEERDQKPASTHSLQTLHTLLSTGSPLKPQSY
EYVYSCIKNNVLLGSISGGTDIISCFMGQNMTVPVYRGEIQARNLGMAVESWSCEGKPVW
GESGELVCLKPIPCQPTHFWNDENGSKYHKAYFSTFPGVWAHGDYCKINPKTGGVVMLGR
SDGTLNPNGVRFGSSEIYNIVEAFDEVSDSLCVPQYNSDGEERVILFLKMGPNKSFSQEL
VGKIRGAIRVALSARHVPALILETKDIPYTISGKKVEVAVKQVIAGKEVTQRGAFSNPDS
LDLYKNLPELQNF''')

class PfamSeqSearch(threading.Thread):
	# Note: Sequences should be added to queue in format:
	# (Unique_ID, Sequence)
	seqQueue = Queue.Queue()
	resultQueue = Queue.Queue()
	SeqData = {}
	seqdata_lock = threading.Lock()

	def __init__(self):
		self.threadid = os.urandom(8)
		self.searchurl = 'http://pfam.sanger.ac.uk/search/sequence'
		threading.Thread.__init__(self)

        """
	def run(self):
		headers = {'Expect':''}
		data_values = {'seq':'', 'output':'xml'}
		# Submit all sequences in queue to Pfam for analysis
		while True:
			try : seqid, seq = PfamSeqSearch.seqQueue.get(block = False)
			except Queue.Empty:
				break
			data_values['seq'] = seq
			data = urllib.urlencode(data_values)
			request = urllib2.Request(self.searchurl,
							data, headers)
			responce = urllib2.urlopen(request)
			responce = xml.dom.minidom.parse(responce)
			results = responce.getElementsByTagName('result_url')
			if not results:
				logging.warning('%s input rejected by \
Pfam or unable to parse responce.' % seqid)
				raise PfamServerError(seqid, self.searchurl, self.threadid, 'None')
			else: result_url = results[0].firstChild.data
			PfamSeqSearch.resultQueue.put((seqid, result_url, \
							self.threadid))

		# Retrieve results from Pfam if available
		while True:
			try: seqid, result_url, threadid = PfamSeqSearch.\
						resultQueue.get(block = True,
								timeout = 5.0)
			except Queue.Empty:
				break

			try: request = urllib2.urlopen(result_url)
			except urllib2.HTTPError as Error:
				request = Error

			if request.code == 202:
				RequeueTimer = threading.Timer(2.0,
					PfamSeqSearch.resultQueue.put,
					args=[(seqid, result_url, threadid)])
				RequeueTimer.start()
			elif request.code == 200:
				PfamSeqSearch.seqdata_lock.acquire()
				if seqid in PfamSeqSearch.SeqData.keys():
					logging.info('Sequence identifier %s \
already contained in PfamSeqSearch.SeqData' % seqid)
				PfamSeqSearch.SeqData[seqid] = \
					PfamSeqSearch.parsepfammatchxml(request)
				if not PfamSeqSearch.SeqData[seqid]:
					logging.info('No domains returned \
for %s' % seqid)
				PfamSeqSearch.seqdata_lock.release()

			else:
				raise PfamServerError(seqid, result_url, \
							threadid, request.code)
        """

	@staticmethod
	def parsepfammatchxml(xml_responce):
		xmldoc = xml.dom.minidom.parse(xml_responce)
		Domains = []
		matches = xmldoc.getElementsByTagName('match')
                matches = [match for match in matches if
                        match.attributes['type'].value == 'Pfam-A']
		for match in matches:
			PfamAccession = match.attributes['accession'].value
			for location in match.childNodes:
				if location.nodeName == 'location':
					if 'significant' in \
						location.attributes.keys() and\
					location.attributes\
						['significant'].value != '1':
							continue
					StartIndex = int(location.attributes\
							['start'].value)
					EndIndex = int(location.attributes\
							['end'].value)
					Evalue = float(location.attributes\
							['evalue'].value)
					Domains.append((PfamAccession,
							StartIndex,
							EndIndex,
							Evalue))
		return sorted(Domains, key = lambda x: x[1])

