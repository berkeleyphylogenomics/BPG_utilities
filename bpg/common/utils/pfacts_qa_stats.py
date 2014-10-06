#!/usr/bin/python
# File: pfacts_qa_stats.py
# Author: Grant Shoffner
# Desc: Prints interesting data from the pfacts database for use in QA.

import sys, logging
import threading, Queue
import optparse
import pfacts003.phylofacts.models as pfacts

class FamilyStats(threading.Thread):
    """Look up basic numbers for families."""
    datalock = threading.Lock()
    data = {}    

    def __init__(self):
        threading.Thread.__init__(self)
    
    def run(self):
        FamilyStats.datalock.acquire()
        FamilyStats.data["Total Families: "] = \
            pfacts.Family.objects.count()
        FamilyStats.data["Families with status 'draft': "] = \
            pfacts.Family.objects.filter(status = 'draft').count()
        FamilyStats.data["Families with status 'bad': "] = \
            pfacts.Family.objects.filter(status = 'bad').count()
        FamilyStats.data["Families with other status: "] = \
            pfacts.Family.objects.exclude(status = 'bad').\
                                exclude(status = 'draft').count()
        FamilyStats.data["Global homology families: "] = \
            pfacts.Family.objects.filter(family_type__id = 'G').count()
        FamilyStats.data["Families with canonical NJ tree: "] = \
            pfacts.Family.objects.filter(
                                    canonical_tree__method = 'nj').count()
        FamilyStats.data["Families with canonical ML tree: "] = \
            pfacts.Family.objects.filter(
                                    canonical_tree__method = 'ml').count()
        FamilyStats.data["Draft families with no NJ/ML tree: "] = \
            pfacts.Family.objects.filter(status = 'draft').\
                exclude(canonical_tree__method = 'ml').\
                exclude(canonical_tree__method = 'nj').count()
        FamilyStats.datalock.release()
        return

class SeqHdrStats(threading.Thread):
    """Stats for the SequenceHeader class."""
    datalock = threading.Lock()
    data = {}

    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        SeqHdrStats.datalock.acquire()
        SeqHdrStats.data["Total SeqHdr objects: "] = \
            pfacts.SequenceHeader.objects.count()
        SeqHdrStats.data["SeqHdr with UniProt accession: "] = \
            pfacts.SequenceHeader.objects.\
                filter(uniprot__accession__contains = '').count()
        SeqHdrStats.data["SeqHdr with sequence chars: "] = \
            pfacts.SequenceHeader.objects.filter(
                sequence__chars__contains = '').count()
        SeqHdrStats.datalock.release()
        return

class UniProtStats(threading.Thread):
    """Stats for the UniProt class."""
    datalock = threading.Lock()
    data = {}

    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        UniProtStats.datalock.acquire()
        UniProtStats.data["Number of UniProt objects: "] = \
            pfacts.UniProt.objects.count()
        UniProtStats.data["SwissProt sequences: "] = \
            pfacts.UniProt.objects.filter(in_swissprot_f = True).count()
        UniProtStats.data["Fragment sequences: "] = \
            pfacts.UniProt.objects.filter(is_fragment = True).count()
        UniProtStats.data["Records with taxonomy: "] = \
            pfacts.UniProt.objects.filter(taxon__id__contains = '').count()
        UniProtStats.datalock.release()
        return

def print_results(classes, output_file):
    """For each class, print their 'data' class variable."""   
    for DataClass in classes:
        DataClass.datalock.acquire()
        output_file.write("--- " + DataClass.__name__ + " ---\n")
        for data in DataClass.data:
            output_file.write(data + repr(DataClass.data[data]) + '\n')
        output_file.write('\n')
        DataClass.datalock.release()
    return

def parse_command_line_options():
    Parser = optparse.OptionParser()
    Parser.add_option('-o', '--output',
                dest = 'output_file',
                default = sys.stdout,
                help = "send output to file, instead of stdout")
    return Parser.parse_args()

def ThreadManager(DataClasses):
    Threads = [DataClass() for DataClass in DataClasses]
    for thread in Threads:
        thread.start()
    for thread in Threads:
        thread.join()
    return

def main():
    """Run all the stats and print to standard out."""
    DataClasses = [FamilyStats, SeqHdrStats, UniProtStats]
    CmdLineOps, Args = parse_command_line_options()
    ThreadManager(DataClasses)
    print_results(DataClasses, CmdLineOps.output_file)
    return

def debugging_tests():
    """Junk code for debugging purposes."""
    logging.warning("Running debugging tests...")
    pass

if __name__ == "__main__":
    main()
    #debugging_tests()
