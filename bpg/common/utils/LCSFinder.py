'''
Finds the longest common string (LCS) in a list of input strings. If run as its 
own script, takes in a file with one string per line. Uses the suffix_tree module
provided in /clusterfs/ohana/sandbox/ijoseph/site-packages/lib64/python/ 
Add the above folder to sys.path to use. 
Prints found LCS.

Created on September 22, 2011

@author: Isaac Joseph (ijoseph@berkeley.edu)
'''
try:
    from suffix_tree import GeneralisedSuffixTree
except ImportError:
    "suffix_tree implementation not found. Make sure /clusterfs/ohana/sandbox/ijoseph/site-packages/lib64/python/ is on path" 
    
import re, sys # built-in Python modules

class LCSFinder(object):
    '''
    Class object to find the least common substring of an input list of strings. 
    '''
    
    def __init__(self, stringList, ignoreCase = False):
        '''
        Constructor; takes in list of strings, whether or not to ignore case when
        finding substrings
        '''
        if ignoreCase:            
            originalList = stringList[:] # Deep copy of list 
            stringList = [s.lower() for s in stringList] # Make lowercase
            
        self.lcs = self.__computeLCS(stringList)[0]          
        
        # If we converted the stringList to lowerCase, search for this string in the original list 
        # (still has correct case)
        if ignoreCase:            
            for string in originalList:                
                m = re.match(self.lcs, string, re.I)
                if m is not None: 
                    self.lcs = m.group(0)
                    return
            return
        
    
    def __computeLCS(self, stringList):
        '''
        Returns a one-element list containing the LCS of the input stringList  
        '''
        
        alphabet = self.__getAlphabet(stringList) # get alphabet of (all characters in) stringList
        
        # check if alphabet requires too many characters to create enough terminal characters 
        # for each string in stringList
        if not self.__isComputable(stringList, alphabet):
            strLstLen = len(stringList)            
            return self.__computeLCS(self.__computeLCS(stringList[0:strLstLen/2]) + \
                                     self.__computeLCS(stringList[strLstLen/2:strLstLen]))
            
        
        (stringList, translationDict) = self.__translateCharacters(stringList, alphabet) # translate characters in stringList
        
        # make suffix tree
        stree = GeneralisedSuffixTree(stringList)        
        # get all shared substrings
        sharedSubstrings = []
        for shared in stree.sharedSubstrings():
            for seq, start, stop in shared:
                sharedSubstrings += [stree.sequences[seq][start:stop]]        

        # find the index of the longest shared substring        
        substringLens = [len(substring) for substring in sharedSubstrings]
        if substringLens == []:
            lcs = [""]
            return lcs
        longestSubstringIndex = substringLens.index(max(substringLens))        
        
        lcs = sharedSubstrings[longestSubstringIndex]
        # Back translate
        for (translatedChar, originalChar) in translationDict.iteritems():
            lcs = lcs.replace(translatedChar, originalChar)
        return [lcs]

    
    def __translateCharacters(self, stringList, alphabet):
        '''
        The characters corresponding to ASCII character code 0 through len(stringList) cannot be included
        in any string in stringList. Thus, we find the set of all characters in the stringList (the
        alphabet of the stringlist), and then translate these to characters whose ASCII code is above
        len(stringList) if possible. 
        '''       
        
        # Creates a dictionary to translate characters in stringList to characters outside of needed terminal character range
        forTransDict = dict(zip(list(alphabet), map(chr, range(255, 255 - len(alphabet), -1))))            
        
        for i in range(len(stringList)):
            for char in forTransDict.keys():
                stringList[i] = stringList[i].replace(char, forTransDict[char])
                
        
        translationDictionary = dict((v,k) for k, v in forTransDict.iteritems())        
        
        
        return (stringList, translationDictionary)
    
    def __getAlphabet(self, stringList):
        '''
        Gets the alphabet of the string list (all characters)
        '''        
        alphabet = set()
        for string in stringList:
            for char in string:
                alphabet = set.union(alphabet, char)
        
        return alphabet
    
    def __isComputable(self, stringList, alphabet):
        '''
        Check to see if we have enough characters left over (if there are more characters needed
        we will need to break up the string list)         
        '''
        return  (255 - len(alphabet) > len(stringList))
            
        
        
    def getCommonName(self):
        return self.lcs

def parseInputFile(filePath):    
    f = open(filePath)
    return [line.rstrip() for line in f.readlines()]

def main():
    assert len(sys.argv) in [2,3], "Usage: LCSFinder <inputFile> <ignoreCase>? ,"+ \
    " where ignoreCase is '1' means case should be ignored"
    
    filePath = sys.argv[1]
    if len(sys.argv) == 3:
        ignoreCase = (sys.argv[2] == '1')            
    else:
        ignoreCase = False
    
    lcsf = LCSFinder(parseInputFile(filePath), ignoreCase)
    print lcsf.getCommonName()

if __name__ == '__main__':
    main()
