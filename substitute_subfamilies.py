#! /usr/bin/python

from optparse import OptionParser
import string, re
from matchmaker.shmm_shmm_lib import *

def subfam_dir(seed_id):
    return os.path.join(cur_dir(), seed_id)

def kerf_dir(seed_id, cutoff):
    return os.path.join(subfam_dir(seed_id), "kerf%d" % cutoff)

def nj_tree_file_path(seed_id):
    return os.path.join(subfam_dir(seed_id), "%s.nj" % seed_id)

def kerf_subfam_file_path(seed_id, cutoff):
    return os.path.join(kerf_dir(seed_id, cutoff), "kerf.trees")

def results_file_path(seed_id, cutoff):
    return os.path.join(subfam_dir(seed_id),
                        "%s_%s_reduced.nj" % (seed_id, cutoff))

def make_reduced_tree(seed_id, cutoff):
    # grab the neighbor joining tree and open the local copy
    nj_tree_file = open("%s" % nj_tree_file_path(seed_id), 'r')
    # open the kerf subfamily file
    subtree_file = open("%s" % kerf_subfam_file_path(seed_id, cutoff), 'r')
    # open the output file 
    supertree_reduced_file = open("%s" % results_file_path(seed_id, cutoff), 'w')
    
    subfam = subtree_file.readline().lstrip("%").rstrip()
    subtree = subtree_file.readline().rstrip().rstrip(";")
    
    line = nj_tree_file.readline()
    supertree = ""
    while (line != ""):
        #line.rstrip()
        line = line[:-1]  # A hack to remove the newline, because rstrip fails
        supertree += line
        line = nj_tree_file.readline()
    
    nj_tree_file.close()
    
    split_pattern = re.compile(":\d")
    insert = ":([\d\.]*)"
    
    while (subfam != ""):
        if (subtree[-1] == ")"):
            subtree = subtree[1:-1]
        subtree_escaped = re.escape(subtree)
        subtree_escaped_split = re.split(split_pattern, subtree_escaped)
    
        subtree_pattern = subtree_escaped_split[0]
        for i in range(1, len(subtree_escaped_split)):
            subtree_pattern = subtree_pattern + insert + subtree_escaped_split[i]
    
        subtree_pattern = re.compile(subtree_pattern)
        m = re.search(subtree_pattern, supertree)
        if (m == None):
            print "WARNING: Subtree not found in original tree."
            print subtree
        else:
            supertree = re.sub(subtree_pattern,
                               subfam+":"+m.group(len(subtree_escaped_split)-1),
                               supertree)
    
        subfam = subtree_file.readline().lstrip("%").rstrip()
        subtree = subtree_file.readline().rstrip().rstrip(";")
    
    subtree_file.close()
    supertree_reduced_file.write(supertree)
    supertree_reduced_file.close()

def main():
    # parse command line options
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--sequence", dest="sequence", 
               help="pdbid for the sequence to build a tree for")
    opt_parser.add_option("-p", "--percentid", type="int", dest="percentid", 
               default=20, help="Minimum %id at which subfamilies were cut")
    (options, args) = opt_parser.parse_args()
    
    if options.sequence is "":
        print "You must indicate a sequence.  Exiting."
        return
    seed_id = options.sequence
    kerf_id = options.percentid
    
    # make sure inputs are valid
    error = False
    if not os.path.exists(subfam_dir(seed_id)):
        print "Cannot find subfamilies for sequence "+seed_id
        error = True
    elif not os.path.exists(kerf_dir(seed_id, kerf_id)):
        print "Cannot find kerf cuts at "+kerf_id+" for sequence "+seed_id
        error = True
    if error:
        print "Aborting due to error.\n"+ \
              "Please ensure the necessary files are precomputed."
        return
    
    print "Pruning subfamilies for sequence %s..." % (seed_id)
    
    # make reduced neighbor joining trees for each seed id
    make_reduced_tree(seed_id, kerf_id)
    
    print "Success.  Output file located in:"
    print "%s" % results_file_path(seed_id, kerf_id)

if __name__ == "__main__":
  main()
