import getopt
import math
import marshal
import operator
import os
import sam as Hmm
import sys
import tempfile
import types

import BPG_common.gnuplot
from bpg.common import fasta
from bpg.common.legacy_utils import *

class dcheuristic:
    def __init__(self, name, v=0, vv=0):
        self.Verbose = v
        self.VeryVerbose = vv
        self.inputmodel = None
        self.mode = 0
        self.threshold = None
        self.generategraph = None
        self.quickview = None
        self.croppedseqs = None
        self.uniqueness = .9
        self.seed = None
        self.name = name
        self.values = None
        self.acceptedoptions = ("u:c:s:hvt:w:g:q", ("help", "user_manual"))
        self.op = operator.and_
        self.doc = "Not set"
        self.isbinary = None

    # seed: file or dictionary type (depending on function that uses it)
    def setSeed(self, seed):
        self.seed = seed

    # model: modelfile to read from (opt)
    def setModel(self, model):
        self.inputmodel = model

    # align: string of alingment file to read
    def setAlignment(self, align):
        self.alignmentfile = align

    # window: window length to look at
    def setWindow(self, window):
        self.window = window

    # threshold: threshold value to set
    def setThreshold(self, threshold):
        self.threshold = threshold
        
    # Input: afile: filename of alignment file to read
    #        seed: file to read seed sequence from
    #        uniqueness: % identity of seqs to remove before generating
    #                    statistics
    # Output: returns a list of floating point numbers of an average columnar
    #         heuristic in the alignment
    def getValues(self, afile=None, seed=None, model=None, uniqueness=None):
        if afile:
            self.setAlignment(afile)
        if seed:
            self.setSeed(seed)
        if model:
            self.setModel(model)
        if uniqueness == None:
            uniqueness = self.uniqueness
        if not self.values:
            self.values = self.generateValues(uniqueness)
        return self.values
    
    # generateValues() should be overridden!
    def generateValues(self, uniqueness):
        if afile:
            self.setAlignment(self.afile)
        if seed:
            self.setSeed(self.seed)
        if model:
            self.setModel(self.inputmodel)
        return None

    def writeSeed(self, filename=None):
        if filename == None:
            filename = "%s.croppedseeds" % self.alignmentfile
        fasta.WriteSequences(seed, filename)

    def main(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], self.acceptedoptions[0], self.acceptedoptions[1])
        except getopt.GetoptError:
            self.usage()
            sys.exit(1)
        for o, a in opts:
            if o in ("-h", "--help"):
                self.usage()
                sys.exit(0)
            elif o == "--user_manual":
                self.user_manual()
                sys.exit(0)
            elif o == "-i":
                self.inputmodel = a
                if not os.access(self.inputmodel, os.F_OK):
                    print "Input Model doesn't exist"
                    sys.exit(3)
            elif o == "-a":
                self.isbinary = 1
            elif o == "-v":
                self.VeryVerbose = 1
                self.Verbose = 1
            elif o == "-w":
                self.window = int(a)
                self.mode = 2
            elif o == "-s":
                self.seed = a            
            elif o == "-t":
                self.threshold = float(a)
                self.mode = 2
            elif o == "-u":
                self.uniqueness = float(a)                
            elif o == "-g":
                self.generategraph = a
            elif o == "-q":
                self.quickview = 1
            elif o == "-c":
                self.croppedseqs = a                
        if len(args) != 1:
            self.usage()            
            sys.exit(5)
        elif self.mode == 2 and (not self.window and self.threshold) or (not 'threshold' in dir(self) and self.window):
            print "Must specify both a window and threshold"
            sys.exit(2)
        elif self.croppedseqs and self.mode != 2:
            print "Cannot write out cropped sequences if no window and threshold specified"
            sys.exit(3)            
        elif self.quickview and self.generategraph:
            print "Cannot show a quickview and generate graph together!"
            sys.exit(4)
        elif self.isbinary and not self.inputmodel:
            print "Need a model to convert to ascii"
            sys.exit(5)

        if self.isbinary:
            self.inputmodel = Hmm.hmmconvert(tempfile.mktemp(), self.inputmodel)            
        values = self.getValues(args[0], self.seed, model = self.inputmodel)
        if self.isbinary:
            os.unlink(self.inputmodel)
        if not values:
            print "Seed and alignment too similar"
            sys.exit()
        for pos in range(len(values)):
            print "%10d\t%10f" % (pos+1, values[pos])

        if self.generategraph:
            type = 'png'            
            gnuplot.CreateGraph(self.generategraph, self.values, xlabel="Columns/Residue Positions", ylabel=self.name, title=args[0])
        elif self.quickview:
            gnuplot.CreateGraph(self.generategraph, self.values, xlabel="Columns/Residue Positions", ylabel=self.name, title=args[0], type='x11')        
        if self.mode == 2:
            print self.doesItFailTheThreshold()
            if self.croppedseqs:
                self.findBestWindowsToRemove(smooth=1)
                sequences, seed = self.cropSelectedHMMRegions(self.croppedseqs)
##                 fasta.WriteSequences(sequences, self.croppedseqs)
##                 self.writeSeed("%s.seed" % self.croppedseqs)

    staticmethod(main) # makes main a static method

    # Print out usage information
    def usage(self):
        print "This will generate statstics and print them to standard out.  Optionally, if both -w and -t are set, then it will return a 1 or 0 depending on whether there is or isn't domain creep."
        print "Usage: %s [options] alignment_file" % sys.argv[0]
        print "Options:"
        print "  -u:\t\t\t%% identity to remove from alignment files (default %.2f)" % self.uniqueness
        print "  -g <outputfile>:\tpng to write"
        print "  -q:\t\t\tquick view - will quickly show the graph"    
        print "  -w:\t\t\twindow size to check for domain creep at any one time"
        print "  -t:\t\t\tthreshold - minimum to not be considered domain creep"
        print "  -s:\t\t\tseed file - file contianing aligned seed"
        print "  -g:\t\t\twill create a graph in png format"
        print "  -c <outputfile>:\tif given window & threshold, will write out new cropped seqs"
        print "  -v:\t\t\tverbose"
        print "  -h,--help:\t\tthis help message"
        print "  --user_manual:\tprint algorithm specifics"

    def user_manual(self):
        print self.doc

    # input: values - a list of averages ordered by position
    #        seqs - hash by id of sequences
    #        window - int of window size to look at
    #        threshold - float of threshold to be <= to be considered domain creep
    #        op - operator fucntion to use
    def doesItFailTheThreshold(self):
        """Look over a given set of averages and check to see if the average in a window of length k is over the threshold.  Return true or false depending on whether it passes the threshold.
        Algorithm:
        - Over a window of length k, take an average of the values
        - compare against threshold, if op threshold, stop, report true
        - else report false"""
        return operator.truth(self.findBestWindowsToRemove())

    # Will check, using a very crude heuristic, where domain creep exists
    def findBestWindowsToRemove(self, smooth=1):
        values = [None]*len(self.values)
        if smooth:
            for x in range(len(self.values)):
                if x < 3 or x > len(self.values)-3:
                    values[x] = self.values[x]
                    continue
                values[x] = (reduce(operator.add, self.values[x-2:x+3]))/float(5)
        else:
            values = self.values
            
        begin = 0
        end = self.window
        assert len(values) >= self.window
        indices = []
        add = operator.add
        while end < len(values) and begin != end:
            average = float(reduce(add, values[begin:end]))/self.window
            if self.op(average, self.threshold):
                begincutback = 0
                endcutback = 0
                # first check the end values to see if they need to be cut or not
                if not self.op(values[begin], self.threshold):
                    while (begin < len(values)) and not self.op(values[begin], self.threshold):
                        begin += 1
                    if not begin < len(values) or not self.op(values[begin], self.threshold):
                        begin -= 1
                    begincutback = 1
                if not self.op(values[end-1], self.threshold):
                    while end >= 0 and not self.op(values[end-1], self.threshold):
                        end -= 1
                    if not end >= 0 or not self.op(values[end-1], self.threshold):
                        end += 1                    
                    endcutback = 1

                # if the ends weren't cut, then try extending
                if not begincutback:
                    while begin-1 >= 0 and (self.op(values[begin-1], self.threshold) or (begin-2 >= 0 and self.op((values[begin-2]+values[begin-1])/float(2), self.threshold)) or (begin-3 >= 0 and self.op((values[begin-3]+values[begin-2]+values[begin-1])/float(3), self.threshold)) or (begin-4 >= 0 and self.op((values[begin-4]+values[begin-3]+values[begin-2]+values[begin-1])/float(4), self.threshold))):
                        begin -= 1
                        
                if not endcutback:
                    while end < len(values) and (self.op(values[end-1], self.threshold) or self.op((values[end-1]+values[end])/float(2), self.threshold)):
                        end += 1
                    if not self.op(values[end-1], self.threshold):
                        end -= 1

                # Now append the range to a list of indices to remove

                if not (end - begin -1) < .7 * self.window:
                    if end > len(values) - self.window:
                        indices.append((begin,len(values)))
                    elif begin < self.window*1.5:
                        indices.append((0, end))
                    elif begin<1.2*self.window or end > len(values)-1.2*self.window:
                        indices.append((begin,end))
                    elif end > len(values)-.2*self.window:
                        indices.append((begin,len(values)))
                        break


                # set the begin and end variables appropriately
                if begincutback:
                    begin += 1
                    end = begin+self.window
                else:
                    begin = end
                    end += self.window

            else:
                begin += 1
                end += 1

            if end > 1.5*self.window and begin < 1.5*self.window:
                begin = max(begin,int(len(values) - 1.5*self.window - 1))
                end = max(end,len(values) - self.window)
            
        # Do some analysis on the self.windows (ie, if within say, 20% of self.window length
        # between each other, join (must be at least 1))
        delta = int(.2 * self.window)
        betterlist = []
        if not delta: delta = 1
        prevpair = [0,0]
        for x,y in indices:
            if x-prevpair[1] <= delta:
                if tuple(prevpair) in betterlist:  betterlist.remove(tuple(prevpair))
                prevpair[1] = y
                betterlist.append(tuple(prevpair))
            else:
                betterlist.append((x,y))
                prevpair = list((x,y))
        self.bestwindows = tuple(betterlist)
        return self.bestwindows

    # alias the name, but do this so we can use it in base classes too
    def isThereDomainCreep(self):
        return self.doesItFailTheThreshold()

    # Crop the selected portions out of the match states of the sequence
    # Input: writeout - write out an aligned file in fasta format of new
    #                   cropped sequences if given a name
    # Output: sequences - hash of sequences by id
    #         seed - cropped seed(s)
    def cropSelectedHMMRegions(self, writeout=0):
        if not self.bestwindows: self.findBestWindowsToRemove()
        seqs = fasta.ReadSequences(self.alignmentfile)
        if isinstance(self.seed, type('')):
            seed = fasta.ReadSequences(self.seed)
        else:
            seed = self.seed
        assert isinstance(seed, type({}))
        for x in seed:
            seed[x] = fasta.matchupperdash.findall(seed[x])
        selectedregions = map(fasta.matchupperdash.findall, seqs.values())
        keys = seqs.keys()
        for x in range(len(keys)):
            selectedregions[x] = [keys[x], selectedregions[x]]
        setslice = operator.setslice
        for x,y in self.bestwindows:
            assert len (seed[seed.keys()[0]]) == len (selectedregions[0][1])        
            Print("Removing regions [%d:%d]" % (x+1,y))
            for j in seed:
                setslice(seed[j], x, y, [('') for p in range(x,y)])
            for seqnum in range(len(selectedregions)):
                setslice(selectedregions[seqnum][1], x, y, [('') for p in range(x,y)])

        sequences = {}
        for x in range(len(selectedregions)):
            sequences[selectedregions[x][0]] = "".join(selectedregions[x][1])
        for x in seed:
            seed[x] = "".join(seed[x])
        if writeout:
            fasta.WriteSequences(sequences, writeout)
            fasta.WriteSequences(seed, "%s.seeds" % writeout)
        return (sequences, seed)

    # This function will take in a tuple of tuple pairs and calculate their
    # true length (since they're in the form of python slices)
    # Output: integer of the length
    def lengthOfRegions(self):
        total = 0
        for x,y in self.bestwindows:
            total += y-x+1
        return total

    # Looks for seed in a couple places (either given, or first seq in
    # alignment) and sets the seed(s) of the object
    # Input: alignmentfile - alignmentfile to possibly look through
    #        matchingregion - list of seqs of matchstates/gaps
    #        seeds - seed object to possibly look through
    #        seqs - a sequence hash by id, previously uniqued
    # Ouput: a tuple(seed id of first seed, seedseq of first id,
    #        matchingregion, seed)
    def findSeed(self, matchingregion=[], seqs={}):
        seeds = {}
        if not self.seed:         # in this case read in first sequence
            seedseq = fasta.ReadOneSequence(self.alignmentfile)
            seqs2 = fasta.ReadSequences(self.alignmentfile)
            seedid = seqs2.keys()[seqs2.values().index(seedseq)]
            seedseq = "".join(fasta.matchupperdash.findall(seedseq))        
            if seedseq in seqs.values():
                del matchingregion[matchingregion.index(seedseq)]
            seeds = {seedid: seedseq}
        else:
            if isinstance(self.seed, types.StringType):
                seeds = fasta.ReadSequences(self.seed)
            else:
                seeds = self.seed
            assert len(seeds)
            seedseq = "".join(fasta.matchupperdash.findall(seeds.values()[0]))
            seedid = seeds.keys()[0]
            if seedseq in matchingregion:
                matchingregion.remove(seedseq)
        assert seeds
        return (seedid, seedseq, matchingregion, seeds)
