'''
This module runs INTREPID but first checks the sequence of interest has been found in the MSA. Note that the actual INTREPID algorithm is implemented in a Perl script (see call below).

It is based heavily on run_intrepid_against_manuscripts.py, found in ohana_repository/intrepid.

Inputs: accession id, path to msa file, path to tree file, path to map file 
Outputs: computed INTREPID scores in "output.aux" file and other output produced by INTREPID.

Author: Curt Hansen
Created: May 14, 2012
Modified: 
'''

import os, sys, glob, pickle, re, shutil
from bpg.common.utils.dir_of_family import get_dir_of_family_accession as get_dir #Need this for Intrepid


def get_file(family_dir,extension):
    '''get the msa and tree file'''
    file=glob.glob(os.path.join(family_dir,extension))
    if len(file) != 1:
        return 'error'
    else:
        return file[0]


def check_inputs(protID,famID):
    OKstatus = [True,True,True] #Initialize to all OK.
    data = []
    family_dir = get_dir(famID)
    msa_file = get_file(family_dir, '*mafft.afa')
    if msa_file=='error':
        OKstatus[0] = False
    else:
        data.append(msa_file)
    tree_file = get_file(family_dir, '*ml.rooted.tre')
    if tree_file=='error':
        OKstatus[1] = False
    else:
        data.append(tree_file)
    idmap_file = get_file(family_dir, '*ungapped.idmap')
    if idmap_file=='error':
        OKstatus[2] = False
    else:
        data.append(idmap_file)
    return [OKstatus,data]


def score_protein(protID,dir,data):
    intrepid_input_name = 'reformatted_MSA.afa'
    MSAfile = data[0]
    treefile = data[1]
    idmapfile = data[2]
    
    h = open(MSAfile)     # read the original MSA
    msa = h.read()
    h.close()

    h = open(idmapfile)   # read the mapping of SEQ#s to headers
    idmap = pickle.load(h)
    h.close()

    seqidFound = False
    for seqno,header in idmap.items():
        msa = msa.replace('>%s\n' % header, '>%s\n' % seqno) #replace the headers with SEQ#s
        if re.search(protID,header):
            seqidFound = True
            seqID = seqno  #record the seqID for the inputID and  
    if not seqidFound:
        return False

    #create dict of header: sequence pairs
    inseqs = dict([i.split('\n',1) for i in msa.strip().lstrip('>').split('\n>')])
    for key,val in inseqs.items():
        inseqs[key] = val.replace('\n', '')     #remove newline chars

    seed_seq = inseqs[seqID] # load the example sequence and remove columns where the seed is gapped
    keys, seqs = zip(*inseqs.items())
    newcols = []
    for num,col in enumerate(zip(*seqs)):
        if not (seed_seq[num] == '-' or seed_seq[num] == '.'):
            newcols.append(col)
    newseqs = map(lambda x: ''.join(x), zip(*newcols))

    # make the directory
    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)

    # write the new alignment file
    h = open(intrepid_input_name, 'w')
    for head,seq in sorted(zip(keys, newseqs)):
        h.write('>%s\n%s\n' % (head, seq))
    h.close()

    # write the config file
    h = open('config.txt', 'w')
    h.write('msa_file %s\ntree_file %s\nsequence_id %s\n' %(intrepid_input_name,treefile,seqID))
    h.close()

    # run intrepid
    os.system('intrepid.pl config.txt 1> intrepid.out 2> intrepid.err')
    os.chdir('..')

    return True


def run_intrepid(protID,famID):

    reldata = check_inputs(protID,"bpg0"+str(famID))
    msaok = reldata[0][0]
    treeok = reldata[0][1]
    idmapok = reldata[0][2]
    scoresok = False #Initialize to false.
    path = './'+protID+'_'+str(famID)+'_intrepid/'
    filepath = path+'output.aux'

    if reldata[0]==[True,True,True]: #Means inputs found for protein, so run INTREPID.
        resultOK = score_protein(protID,path,reldata[1]) #Run INTREPID for protein/family combo.
        if resultOK: #Means INTREPID completed successfully.
            try:
                score_file = open(filepath,'r')
                scoresok = True
            except:
                scoresok = False

        if scoresok:
            readFirstLine = False
            scores = []
            for line in score_file:
                if not readFirstLine:
                    readFirstLine = True #Skip the header row.
                else:
                    details = line.rstrip('\n') #Drop trailing carriage return for last field.
                    details = line.split('|') #Parse line.
                    scores.append(details) #Append scores to list of scores.
            score_file.close()
            shutil.rmtree(path,ignore_errors=True) #Delete directory of Intrepid results for protein.

    if scoresok: # At this point, if scoresok is true then INTREPID run is OK.
        return [True, scores]
    else:
        return [False, [reldata[0][0],reldata[0][1],reldata[0][2],scoresok]]
