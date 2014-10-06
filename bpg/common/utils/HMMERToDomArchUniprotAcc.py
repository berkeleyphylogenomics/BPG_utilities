#!/usr/bin/env python2.6
'''
Converts hmmscan --domtblout output from  multiple UniProt sequences into DomainArchitecture-> UniProt Accession table 

Note: requires PFAM-A architecture Pfam scan domain table output as input. For example, the 
file path $2 in the below command ($1 is the fasta file on which hmmscan is run). 

hmmscan --cut_ga --noali --acc --domtblout $2 /clusterfs/ohana/external/pfam/current/Pfam-A.hmm $1 >/dev/null

Created on Dec 7, 2011

@author: Isaac Joseph (ijoseph@berkeley.edu)
'''
import re, sys

class Domain(object):
    '''
    Holds information about a domain: name, accession, description, alignment start, alignment end 
    on specific architecture to which it is associated   
    ''' 
    
    def __init__(self, pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd):
        self.pfamDomainName = pfamDomainName
        self.pfamDomainAcc = pfamDomainAcc
        self.pfamDomainDesc = pfamDomainDesc
        self.aliStart = aliStart
        self.aliEnd = aliEnd
        
    def __eq__(self, other):
        '''
        Override equals method so that we can check if a domain has been
        inserted into an architecture before 
        '''
        if isinstance(other, Domain):
            return self.pfamDomainName == other.pfamDomainName \
                and self.pfamDomainAcc == other.pfamDomainAcc \
                and self.pfamDomainDesc == other.pfamDomainDesc \
                and self.aliStart == other.aliStart \
                and self.aliEnd == other.aliEnd
        else: 
            return False
    
    def getVersionlessAccession(self):
        versionlessAcccession = re.search('PF[^\.]+', self.pfamDomainAcc).group(0)
        if versionlessAcccession:
            return versionlessAcccession
        else:
            return self.pfamDomainAcc

class DomainArchitecture(object):
    '''
    Holds a protein accession's domain architecture, which consists of a collection of 
    Domains
    '''

    def __init__(self, uniprotAccession):
        '''
        Takes in a protein accession. 
        '''
        self.uniprotAccession = uniprotAccession        
        self.domainList = [] 
        self.numOverlappingDomains = 0 # number of domains that overlap
        
    def addDomain(self, pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd):
        '''
        Adds a domain to the architecture for this protein if another domain is not already extant 
        in the same exact location.
        '''
        newDomain = Domain(pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd)
        if newDomain not in self.domainList:
            self.domainList += [newDomain]        
    
    def retrieveOrderedDomainAccessions(self):
        ''' 
        Orders all Domain objects in list by first alignment index; then checks
        for domains that have the start alignment position, and orders by end alignment.        
        '''        
        def orderDomainCompareFn(domain1, domain2):
            ''' 
            Compare the ordering of the two domains based on where they were 
            aligned 
            '''            
            if domain1.aliStart != domain2.aliStart: # If start in different places, take first
                return domain1.aliStart - domain2.aliStart
            else:
                self.numOverlappingDomains +=1 # we have an overlapping sequence because they start in the same place 
                return domain1.aliEnd - domain2.aliEnd
        
        self.domainList.sort(cmp = orderDomainCompareFn)
        
        return [d.getVersionlessAccession() for d in self.domainList]
    
    def __str__(self, *args, **kwargs):
        ''' 
        Outputs dash-delimited string of domain architecture
        '''
        return "-".join(self.retrieveOrderedDomainAccessions())    
        
def parseHMMERFile(domainTableFile):
    '''
    Returns hashmap <String uniProtAccession> -> <DomainArchitecture domainArchitecture> 
    '''
    # Number of overlapping domains    
    uniprotAccToDomArch = {}
    
    domTblF = open(sys.argv[1]) # open domain table file
    for line in domTblF:
        if line.startswith('#'): # ignore comment line
            continue
                        
        lineAr = re.split(r'\s+', line) # gather data in line 
        (pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd) = \
        (lineAr[0], lineAr[1], " ".join(lineAr[22:len(lineAr)]).rstrip(), int(lineAr[17]), int(lineAr[18]))
        
        uniprotAccession = lineAr[3]; # grab UniProt accession from this
        uniprotAccession = re.search(r'\|([^\|]+)\|', uniprotAccession).group(1)
        
        # Add the domain associated with this line to the domain architecture for this UniProt accession        
        if uniprotAccession not in uniprotAccToDomArch.keys():
            uniprotAccToDomArch[uniprotAccession] =  DomainArchitecture(uniprotAccession)
                    
        uniprotAccToDomArch[uniprotAccession].addDomain(pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd)
    domTblF.close()
    return uniprotAccToDomArch

def reverseUniProtToArch(inHash):
    domArchToUniProt = {}
    # Reverse hash to go from DomainArchitecutre -> List of UniprotAcc
    for (UniProtAcc, domArch) in inHash.iteritems():
        try:
            domArchToUniProt[str(domArch)] += [UniProtAcc]        
        except KeyError:
            domArchToUniProt[str(domArch)] = [UniProtAcc]    
    return domArchToUniProt

if __name__ == '__main__':
    assert len(sys.argv) == 2, 'Usage: HMMERParser <hmmer domain table output file' + \
    '(hmmscan with --domtblout option)>'    
    # Get domain architecture Strings for each protein    
    uniProtToDomArch = parseHMMERFile(sys.argv[1]) # Gets UniprotAcc -> DomainArchitecture Hash 
    domArchToUniProt = reverseUniProtToArch(uniProtToDomArch)
    
    for domArch in domArchToUniProt.keys():
        for uniProtAcc in domArchToUniProt.get(domArch):
            print domArch, "\t", uniProtAcc
        
    

    
    
