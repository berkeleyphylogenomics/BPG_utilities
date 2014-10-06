'''
This is a module of commands used to process HMMS for Daphne (and any other programs needing them).
Author: Curt Hansen
Created: May 22, 2012
Modified:
'''


import subprocess, os


def hmmbuild(msa_file_name, name):
    args = ["/clusterfs/ohana/software/bin/hmm3build", 
        "-o", "/dev/null",
        "--informat", "afa",
        "-n", str(name),
        "/dev/stdout", 
        str(msa_file_name)]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate()[0]


def hmmsearch(hmm_model, eval):
    args = ["/clusterfs/ohana/software/bin/hmm3search",
            "-E", str(eval),
            "-o", "/dev/null",
            "--tblout", "/dev/stdout",
            "-", "/clusterfs/ohana/external/UniProt/current/protein"]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(hmm_model)[0]


def getfasta(accessions): 
    args = ["/clusterfs/ohana/software/bin/fastacmd",
            "-i", "/dev/stdin",
            "-d", "/clusterfs/ohana/external/UniProt/current/protein",
            "-o", "/dev/stdout"]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(accessions)[0]


def create_fasta_file(accessions,output_file_name,tempdir):
    if not os.path.exists(tempdir):
        os.mkdir(tempdir)
    os.chdir(tempdir)
    args = ["/clusterfs/ohana/software/bin/fastacmd",
            "-i", "/dev/stdin",
            "-d", "/clusterfs/ohana/external/UniProt/current/protein",
            "-o", str(output_file_name)]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,sterr = process.communicate(accessions)
    os.chdir('..')
    

def hmmalign(hmm_model,fasta_file_name,tempdir):
    target = str(tempdir)+'/'+str(fasta_file_name)
    args = ["/clusterfs/ohana/software/bin/hmm3align",
            "-", str(target)]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(hmm_model)[0]

