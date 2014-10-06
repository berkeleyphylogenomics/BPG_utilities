#!/clusterfs/ohana/software/bin/python2.6 -Wignore::DeprecationWarning
'''
This script runs the Daphne algorithm on a pre-computed MSA (created using, for example, MAFFT).

Inputs: "msa" file in afa format
Outputs: Set of scores, one for each cluster of sequences in the family.

Author: Curt Hansen
Created: June 15, 2012
Modified: 
'''

import hmm_functions as h
import hcluster as hc
from Bio import AlignIO
from Bio import Cluster
from cStringIO import StringIO
from datetime import datetime
import sys, os, re
import numpy as n


def main():
    if len(sys.argv) != 4:
        print "\nProgram terminated...!\nUsage: %s <msa_file> <num_clusters> <eval> \n" % sys.argv[0]
        sys.exit(1)

    msa_file = sys.argv[1]
    cntClusters = int(sys.argv[2])
    eval = sys.argv[3]

    if not os.path.exists(msa_file):
        print "\nError: Cannot find file '%s'\nProgram terminated...!\n" % msa_file
        sys.exit(1)
        
    l1 = 65
    l2 = 15

    AAuppercase = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    AAlowercase = ['a','r','n','d','c','e','q','g','h','i','l','k','m','f','p','s','t','w','y','v']
    AAgaps = ['-','.','x','X']
    AAlist = AAuppercase + AAlowercase + AAgaps

    startTime = datetime.now()
    pgmName = 'run_daphne_on_family'
    fout = open(pgmName+".out."+startTime.strftime("%Y%m%d.%Hh%M"),"w")
    fout.write("\nStarting Daphne program...\n")
    fout.write("Turned off deprecation warnings.\n")
    fout.write("\nSteps".ljust(l1)+"Elapsed Time".rjust(l2)+"\n")
    fout.write("-"*(l1+l2)+"\n")

    '''Build HMM '''
    ncolsMSA = AlignIO.read(open(msa_file,'r'),'fasta').get_alignment_length()
    lastTime = startTime
    hmm = h.hmmbuild(msa_file,'hmm_'+msa_file)
    r = re.compile(r'LENG\s+(\d*)\n')
    ncolsHMM = r.findall(hmm)[0]
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Built HMM...".ljust(l1)+str(diff)+"\n")

    '''Search sequence db against HMM'''
    lastTime = nowTime
    hits = h.hmmsearch(hmm,eval)
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Searched sequence database against HMM...".ljust(l1)+str(diff)+"\n")

    '''Obtain FASTA sequences for hits'''
    lastTime = nowTime
    p = re.compile(r"tr\|(\w{6,6})\|")
    accessions = p.findall(hits)
    fasta_file_name = msa_file+'.hmmsearch.fasta'
    tempdir = 'temp.dir.'+msa_file
    h.create_fasta_file('\n'.join(accessions)+'\n',fasta_file_name,tempdir)
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Obtained FASTA sequences for hits in HMM search...".ljust(l1)+str(diff)+"\n")

    '''Align sequences against HMM'''
    lastTime = nowTime
    hmm_alignment_temp = StringIO(h.hmmalign(hmm,fasta_file_name,tempdir))
    hmm_alignment = AlignIO.read(hmm_alignment_temp,'stockholm')
    nrowshmma = len(hmm_alignment)
    ncolshmma = hmm_alignment.get_alignment_length()
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Aligned selected sequences against HMM...".ljust(l1)+str(diff)+"\n")
    
    '''Extract data from alignment'''
    lastTime = nowTime
    list_alignment = []
    for record in hmm_alignment:
        list_alignment.append(list(record.seq))
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Extracted data from alignment...".ljust(l1)+str(diff)+"\n")

    '''Compute scores per entry in alignment'''
    lastTime = nowTime
    counts = n.zeros((len(AAlist),ncolshmma))
    for i in range(nrowshmma):
        for j in range(ncolshmma):
            counts[AAlist.index(list_alignment[i][j]),j] += 1 #numpy uses 0 based indexing
    percents = counts/sum(counts)
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Computed scores per entry in alignment...".ljust(l1)+str(diff)+"\n")
    
    '''Compute clusters and average cluster score profile and dimenstions'''
    lastTime = nowTime
    tree = hc.ClusterTree(percents)
    groups = tree.get_cluster_groups(cntClusters)
    group_aves = n.zeros((cntClusters,ncolshmma))
    for i in range(cntClusters):
        group_aves[i,:] = n.average(percents[groups[i]],0)
    n.savetxt('daphne_profile_scores.csv',group_aves,fmt='%5.4f',delimiter=',')
    nowTime = datetime.now()
    diff = nowTime - lastTime
    fout.write("Computed clusters and average cluster score profiles...".ljust(l1)+str(diff)+"\n")

    fout.write("\nResults\n")
    fout.write("-"*(l1+l2)+"\n")
    fout.write("Number of clusters created:".ljust(l1)+str(cntClusters).rjust(l2)+"\n")
    fout.write("Number of columns in original MSA:".ljust(l1)+str(ncolsMSA).rjust(l2)+"\n")
    fout.write("Number of columns in HMM model:".ljust(l1)+str(ncolsHMM).rjust(l2)+"\n")
    fout.write("Number of columns in HMM MSA:".ljust(l1)+str(ncolshmma).rjust(l2)+"\n")
    endTime=datetime.now()
    diff=endTime-startTime
    fout.write("Elapsed time (HH:MM:SS):".ljust(l1)+str(diff).rjust(l2)+"\n\n")
    fout.close()


if __name__=='__main__':
  main()
