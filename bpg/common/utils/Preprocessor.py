'''
Preprocesses input file of protein names (one per line).

Prints out the preprocessed form of each protein name (see individual classes
in this file for specific details of what is removed). 

Created on Oct 5, 2011

@author: Isaac Joseph (ijoseph@berkeley.edu)
'''
import re, sys
class Preprocessor(object):
    '''
    Superclass for all classes that change a list of words into 
    a modified list of words 
    '''
    
    def __init__(self, nameList):
        self.nameList = nameList   
    
    def getModifiedList(self):
        self.preprocess()
        return self.nameList

    def preprocess(self):
        pass
    
class SymbolRemover(Preprocessor):
    '''
    Removes specified symbols from names: 
    specifically, "(", ",", ";", "-". 
    Doesn't remove dashes if they are of the form 
    number-dash-word
    '''
    def __init__(self, wordList, removedSymbolString = '"(),;_|:'):        
        Preprocessor.__init__(self, wordList)        
        self.removedSymbolString = removedSymbolString
        
    def preprocess(self):        
        for char in self.removedSymbolString:
            if char == "_":
                self.nameList = [word.replace(char, " ") for word in self.nameList]
            else:
                self.nameList = [word.replace(char, "") for word in self.nameList]

class LowMeaningWordRemover(Preprocessor):
    ''' 
    Removes specified words from WordList; if word does not exist in WordList afterwards, removes the entire 
    entry.
    
    '''
    def __init__(self, nameList, wordsToRemove = ["predicted", "putative", \
                                                "uncharacterized", "probable", "function"] ):
        Preprocessor.__init__(self, nameList)        
        self.wordsToRemove = wordsToRemove
        
    
    def preprocess(self):
        '''
        Removes specified words, and remove any entries containing the word 
        'uncharacterized' or 'undetermined'
        '''
        for nameIndex in range(len(self.nameList)):
            wordList = self.nameList[nameIndex].split()            
            for word in wordList:
                if word.lower() in self.wordsToRemove:
                    self.nameList[nameIndex] = self.nameList[nameIndex].replace(word,"")
                    # Remove double and trailing spaces resulting from this
                    self.nameList[nameIndex] = self.nameList[nameIndex].strip()
                    self.nameList[nameIndex] = self.nameList[nameIndex].replace("  ", " ")
                    
class LowMeaningPhraseRemover(Preprocessor):                
    '''
    Removes low-meaning phrases, such as transcript variant
    information, anything before the phrase "similar to", and 
    isoform information, and the the end of the phrase "protein"
    '''
    
    def preprocess(self):
        for nameIndex in range(len(self.nameList)):
            # Remove " isoform" and the word directly following
            m1 = re.search(r'\w*\-isoform|,?\sisoform\s[^\s]+', self.nameList[nameIndex])
            if m1:
                self.nameList[nameIndex] = self.nameList[nameIndex].replace(m1.group(0), '')
            
            # If ends with 'protein', kill that
            if self.nameList[nameIndex].endswith("protein"):
                self.nameList[nameIndex] = self.nameList[nameIndex].rstrip("protein")

            m2 = re.search(r'.*similar to\s', self.nameList[nameIndex], re.I)
            if m2:
                self.nameList[nameIndex] = self.nameList[nameIndex].replace(m2.group(0), '')
                
            # remove everything regarding ", partially confirmed by transcript evidence
            m3 = re.search(r',?\spartially confirmed by transcript evidence', self.nameList[nameIndex])
            if m3:
                self.nameList[nameIndex] = self.nameList[nameIndex].replace(m3.group(0),'')

class IncredulousNameRemover(Preprocessor):
    '''
    Remove (replace with empty string) names that include 
    words that indicate the name is invalid, such as "uncharacterized" or "undetermined"
    '''
    def __init__(self, nameList, invalidatingWordList = \
                 ["uncharacterized", "undetermined", "unassigned", "unknown", \
                  "genome", "chromosome", "sequence"]):
        Preprocessor.__init__(self, nameList)
        self.invalidatingWordList = invalidatingWordList
        
    def preprocess(self):
        '''
        Remove names with specified words
        '''
        for nameIndex in range(len(self.nameList)):
            for invalidatingName in self.invalidatingWordList:
                if invalidatingName.lower() in self.nameList[nameIndex].lower():
                    self.nameList[nameIndex] = ''
        
class UniProtRegexRemover(Preprocessor):
    ''' 
    Removes all words from list that match regexes indicating their likely 
    lack of use, accounting for the possibility of "protein" as a suffix. 
    Also remove entire names that contain the word "uncharacterized" 
    '''
    def preprocess(self):
        for nameIndex in range(len(self.nameList)):
            # remove anything that's one word and ends with at least 3 numbers and a 'p', or 
            # 4 uppercase letters or numbers and then lowercase 'p'
            if len(self.nameList[nameIndex].lower().replace("protein",'').split()) == 1: 
                if  re.search(r'[A-Z0-9]{4}p|\d\d\dp$', self.nameList[nameIndex]):
                    self.nameList[nameIndex] = ''
                    continue
                # remove anything that includes at least 4 consecutive digits 
                if re.search(r'[0-9][0-9][0-9][0-9]', self.nameList[nameIndex]):
                    self.nameList[nameIndex] = ''
                    continue
                
                # remove anything that begins with "Zgc" and includes at least 3 numbers
                if re.search(r'zgc:?[0-9][0-9][0-9]', self.nameList[nameIndex], re.IGNORECASE):
                    self.nameList[nameIndex] = ''
                    continue
            # remove anything that contains uppercase "CDS" as is own word 
            if re.search(r'[^a-zA-Z]CDS[^a-zA-Z]|CDS[^a-zA-Z]|[^a-zA-Z]CDS', self.nameList[nameIndex]):
                if len(self.nameList[nameIndex].strip()) > 3:
                    self.nameList[nameIndex] = ''
                    continue
            
def preprocessAll(nameList):    
    '''
    Preprocesses all elements in input nameList
    '''
    # Make nameList a list if it is just a single name    
    if not isinstance(nameList, list): 
        nameList = [nameList]
    return UniProtRegexRemover(\
                               IncredulousNameRemover(\
                                                      LowMeaningPhraseRemover(\
                                                                              LowMeaningWordRemover(\
                                                                                                    SymbolRemover(nameList).getModifiedList()).getModifiedList()).getModifiedList()).getModifiedList()).getModifiedList()
            
def main():            
    f = open(sys.argv[1])    
    for line in preprocessAll([line.rstrip() for line in f.readlines()]):
        print line

if __name__ == '__main__':
    assert len(sys.argv) == 2, "Usage: Preprocessor <inputFile>"
    main()
    
