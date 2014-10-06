'''
Gets longest common domain architecture of all sequences by parsing -domtblout from hmmscan. 

Accepts as input a domain table output file from hmmscan and outputs the longest common 
subsequence of Pfam domains (Domain table output file usually created with all the UniProt sequences in one PHOG)

Note: requires PFAM-A architecture Pfam scan domain table output as input. For example, the 
file path $2 in the below command ($1 is the fasta file on which hmmscan is run). 

hmmscan --cut_ga --noali --acc --domtblout $2 /clusterfs/ohana/external/pfam/current/Pfam-A.hmm $1 >/dev/null

Created on Nov 15, 2011
@author: Isaac Joseph (ijoseph@berkeley.edu)

'''
import sys, re
from LCSFinder import LCSFinder

class DomainArchitecture(object):
    '''
    Holds a protein accession's domain architecture.
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
        Adds a domain to the architecture for this protein.
        '''
        newDomain = Domain(pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd)
        if newDomain not in self.domainList:
            self.domainList += [newDomain]        
    
    def retrieveOrderedDomainAccessions(self):
        ''' 
        Orders all Domain objects in list by first alignment index; then checks
        for domains that have the start alignment position, and orders by end alignment 
        '''        
        def orderDomainCompareFn(domain1, domain2):
            ''' 
            Compare the ordering of the two domains
            '''            
            if domain1.aliStart != domain2.aliStart: # If start in different places, take first
                return domain1.aliStart - domain2.aliStart
            else:
                self.numOverlappingDomains +=1 # we have an overlapping sequence because they start in the same place 
                return domain1.aliEnd - domain2.aliEnd
        
        self.domainList.sort(cmp = orderDomainCompareFn)
        
        return [d.getVersionlessAccession() for d in self.domainList]

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
                    
        uniprotAccToDomArch[uniprotAccession].addDomain\
        (pfamDomainName, pfamDomainAcc, pfamDomainDesc, aliStart, aliEnd)
    domTblF.close()
    return uniprotAccToDomArch

def getDomainArchitectureStrings(domainArchitectures):
    '''
    Given a hashmap <String uniProtAccession> -> <DomainArchitecture domainArchitecture>, 
    output list [String pfamDomainArchitectureString], removing trailing version numbers on 
    Pfam domains if present. pfamDomainArchitectureString is a concatenation of all 
    Pfam domains in the protein's domain architecture in order.  
    '''
    outList = []
    for domArch in domainArchitectures.values():
        outString = ""
        for domArchAcc in domArch.retrieveOrderedDomainAccessions():
            outString += domArchAcc                
        outList += [outString]
    return outList

def stripPfamString(pfamString):
    '''
    For a common domain architecture string, trims the end to 
    remove fragments of Pfam- A domain architectures (assuming form PFxxxxx) 
    e.g. 45PF12345PF1 -> PF12345; allows us to use simple LCS on characters
    '''
    pfamString = re.sub('^\d+', '', pfamString) # kill leading digits
    # kill trailing stuff if doesn't end with 5 digits (isn't well formed)
    if not re.search('\d{5}$', pfamString): 
        pfamString = re.sub('PF\d*$', '', pfamString)    
    return pfamString

def delimitPfamString(pfamString):
    '''
    Once again assuming PFxxxxx format for PFAM domains, delimit result string
    with dashes
    '''
    separatedDomains = re.split(r'P', pfamString)    
    return "-".join(["P" + dom for dom in separatedDomains][1:len(separatedDomains)])

def main():
    assert len(sys.argv) == 2, 'Usage: HMMERParser <hmmer domain table output file>'    
    # Get domain architecture Strings for each protein
    domainArchitectures = parseHMMERFile(sys.argv[1]) 
    domArchStrings = getDomainArchitectureStrings(domainArchitectures) 
    # Find longest common architecture for each string
    longestCommonArchitectureString = LCSFinder(domArchStrings).getCommonName()
    print delimitPfamString(stripPfamString(longestCommonArchitectureString))
    
if __name__ == '__main__':
    main()

                        
    

    
    
