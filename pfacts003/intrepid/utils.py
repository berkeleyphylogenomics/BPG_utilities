''' This file contains utilities used for intrepid. '''
import os, sys, gzip, cStringIO, re
import subprocess
from pfacts003.utils.hmm import hmmbuild, parse_domtblout
from pfacts003.intrepid.consts import *
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.PDB import PDBParser, Selection
from Bio.PDB.Polypeptide import three_to_one
from Bio.Alphabet import SingleLetterAlphabet

def repairPDB(input_pdb, options=[]):
    args = ['/clusterfs/ohana/software/intrepid/scripts/repairPDB', '/dev/stdin'] + options
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(input_pdb)[0]

def hmmsearch(hmm_file, seq_db, domtbl_file, options):
    args = ['/clusterfs/ohana/software/bin/hmm3search',
        '-o', '/dev/null', '--notextw', '--domtblout',
        domtbl_file] + options + [hmm_file, seq_db]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.communicate()
    return None

def hmmbuild(msa, name, options=[]):
    args = ["/clusterfs/ohana/software/bin/hmmbuild", "-o", "/dev/null", \
        "-n", str(name)] + options + ["/dev/stdout", msa]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate()[0]

def intrepid(cfg, out):
    ''' takes the config file path as input, runs intrepid, output stuff is in directory '''
    args = ['/clusterfs/ohana/software/test/repo_bin/intrepid.pl', cfg]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    f = open(out, 'wb')
    o = process.communicate()[0]
    f.write(o)
    f.close()
    return f.name

def fast_tree(alignment):
    args = ["/clusterfs/ohana/software/bin/FastTree", "-quiet", "-nopr", alignment]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate()[0]

def get_query_rec(msa, qid='QUERY'):
    ''' returns the query record from the msa. '''
    for rec in msa:
        if qid in rec.id:
            return rec
    return None

def get_msa_rec(msa, acc):
    ''' returns the correct record from the msa. '''
    for rec in msa:
        try:
            this_acc = rec.id.split('|')[1]
        except:
            this_acc = ''
        if this_acc == acc:
            return rec
    return None

def convert_a2m(ali):
    fh = cStringIO.StringIO(ali)
    msa = AlignIO.read(fh, 'fasta')
    fh.close()
    new_msa = []
    for rec in msa:
        new_seq = Seq(re.sub(r'[a-z.]', '', str(rec.seq)), SingleLetterAlphabet())
        new_rec = rec
        new_rec.seq = new_seq
        new_msa.append(new_rec)
    new_msa = MultipleSeqAlignment(new_msa)
    return new_msa.format('fasta')
        
def hmmalign(fasta, hmm, options=[]):
    args = ['/clusterfs/ohana/software/bin/hmm3align'] + options + [hmm, fasta]

    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return process.communicate()[0]

def mafft_files(fasta, options=[]):
    '''wrapper for mafft'''
    args = ["/clusterfs/ohana/software/bin/mafft"] + options + [fasta]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate()[0]

def mafft(fasta, options=[]):
    '''wrapper for mafft'''
    args = ["/clusterfs/ohana/software/bin/mafft"] + options + ["/dev/stdin"]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate(fasta)[0]

def get_fasta_for_accession(acc, blastable_db = '/clusterfs/ohana/external/UniProt/2013-03-19/protein'):
    ''' wrapper for fastacmd '''
    args = ['/clusterfs/ohana/software/bin/fastacmd', '-d', blastable_db, '-s', acc]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate()[0]

def new_get_fasta_for_accession(acc, blastable_db = '/clusterfs/ohana/external/UniProt/2013-12-17/protein'):
    ''' wrapper for blastdbcmd. '''
    args = ['/clusterfs/ohana/software/bin/blastdbcmd', '-db', blastable_db, '-dbtype', 'prot', '-entry', acc]
    process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.communicate()[0]

def discern(seed, pdbid, path, offset=0, has_upload=False):
    '''wrapper for discern'''
    os.chdir(path)
    if has_upload:
        args = ['/clusterfs/ohana/software/intrepid/scripts/discern.py', '-s', seed, '-c', pdbid, '-o', '%d' % (-1*offset)]
    else:
        args = ['/clusterfs/ohana/software/intrepid/scripts/discern.py', '-s', seed, '-p', pdbid, '-o', '%d' % (-1*offset)]
    s = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    s.communicate()
    return None

def get_pdb_file_for_id(id, library_path):
    ''' returns the path for the pdb file if it exists '''
    if os.path.exists(os.path.join(library_path, '%s/pdb%s.ent.gz' % (id[1:3],id))):
        return os.path.join(library_path, '%s/pdb%s.ent.gz' % (id[1:3],id))
    return None

def read_chains_from_pdb_file(pid, library_path):
    ''' this function will return the list of chains in the first model of the structure '''
    # first get and unzip the pdb file
    f = gzip.open(get_pdb_file_for_id(pid, library_path), 'rb')
    fh = cStringIO.StringIO(f.read())
    f.close()
    # next set parser options
    parser = PDBParser(PERMISSIVE=True, get_header=True)
    # get the structure
    structure = parser.get_structure(pid, fh)
    # close the file handle
    fh.close()
    # return the chains
    return [r.id for r in Selection.unfold_entities(structure.get_list()[0], 'C')]

def get_chains_from_any_pdb_file(pdb_file, pdb_id):
    ''' this function returns the chains from any pdb file. '''
    pdb_f = cStringIO.StringIO(pdb_file)
    # set the parser
    parser = PDBParser(PERMISSIVE=True, get_header=True)
    # get the structure
    structure = parser.get_structure(pdb_id, pdb_f)
    # return the chain ids for model 0 always
    chains = [c.id for c in structure.get_list()[0].get_list()]
    pdb_f.close()
    return chains

def read_sequence_from_any_pdb_file(pdb_file, pdb_id, chainid):
    ''' this function returns the sequence and position from the pdb file. '''
    pdb_f = cStringIO.StringIO(pdb_file)
    # set the parser
    parser = PDBParser(PERMISSIVE=True, get_header=True)
    # get the structure
    structure = parser.get_structure(pdb_id, pdb_f)
    pdb_f.close()    
    chain = structure[0][chainid]
    # get the residues
    res = [r for r in chain.get_list() if not r.id[0] or r.id[0] == ' ']
    # get the amino acid array
    aa_array = [three_to_one(r.resname) for r in res]
    # get the tuple array
    tuple_array = [( r.id[1], three_to_one(r.resname) ) for r in res]
    return (''.join(aa_array), tuple_array)

def read_sequence_from_pdb_file(pid, chainid, library_path):
    ''' this function will read the sequence and position from the pdb file. '''
    # first get and unzip the relevant pdb file
    f = gzip.open(get_pdb_file_for_id(pid, library_path), 'rb')
    old_pdb = f.read()
    # repair the pdb here
    #args = ['/clusterfs/ohana/software/intrepid/scripts/repairPDB', '-', '-rres']
    #process = subprocess.Popen(args, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #fh = cStringIO.StringIO(process.communicate(old_pdb)[0])
    fh = cStringIO.StringIO(old_pdb)
    f.close()
    # next set parser options
    parser = PDBParser(PERMISSIVE=True, get_header=True)
    # next get structure
    structure = parser.get_structure(pid, fh)
    # close the file handle
    fh.close()
    # get the correct chain - user model 0 for now....might change this in the future
    chain = structure[0][chainid]
    # get the residues from the chain
    res = [r for r in chain.get_list() if not r.id[0] or r.id[0] == ' ']
    # get the array for the return string
    aa_array = [three_to_one(r.resname) for r in res]
    # get the array of tuples for return
    tuple_array = [( r.id[1], three_to_one(r.resname) ) for r in res]
    return (''.join(aa_array), tuple_array) 

            
def get_pdb_fasta_for_id(id, chainid, library_path):
    ''' This returns the fasta for the given pdb id.  If chain id is blank, this will 
        look for the (now longest) (commented out: most represented) chain.  It returns a 
        tuple with the fasta entry and the "name" and the number of members.'''
    if chainid:
        f = new_get_fasta_for_accession(id+chainid, os.path.join(library_path, 'blastdbs', 'pdb'))
        if not f.startswith('>'):
            return None
        return ( f, chainid, 1 )
    else:
        chains = read_chains_from_pdb_file(id, library_path) 
        d = {}
        for chain in sorted(chains, key=lambda x: x):
            f = new_get_fasta_for_accession(id+chain, os.path.join(library_path, 'blastdbs', 'pdb'))
            if f.startswith('>'):
                fseq = ''.join(f.split('\n')[1:])
                if fseq in d:
                    d[fseq]['count'] += 1
                    d[fseq]['name'] += ', %s' % (chain)
                else:
                    d[fseq] = {}
                    d[fseq]['count'] = 1
                    d[fseq]['name'] = '%s' % (chain) 
                    d[fseq]['header'] = f.split('\n')[0]
        # Find the one with the most repeats
        #sd = sorted(d.items(), key=lambda x: x[1]['count'], reverse=True)
        #if sd[0][1]['count'] == 1:
        # If they are all equal, take the longest one    
        if not d:
            return None
        sd = sorted(d.items(), key=lambda x: len(x[0]), reverse=True)
        return ( sd[0][1]['header'] + '\n%s' % (sd[0][0]), sd[0][1]['name'], sd[0][1]['count'] )

def add_reference_annotation_line(input_msa_path, output_msa_path):
    # this assumes a fasta
    this_input_file = open(input_msa_path, 'rb')
    this_msa = AlignIO.read(this_input_file, 'fasta')
    this_input_file.close()
    this_query_rec = get_query_rec(this_msa)
    this_query_string = str(this_query_rec.seq)
    this_annotation_line = '#=GC RF '
    for ch in this_query_string:
        if ch.isupper() or ch == '-':
            this_annotation_line += 'x'
        else:
            this_annotation_line += '.'
    this_annotation_line += '\n'
    these_stockholm_lines = this_msa.format('stockholm')
    output_msa = ''
    for line in these_stockholm_lines.split('\n'):
        if line:
            if line.startswith('//'):
                #write the annotation line at the end
                output_msa += this_annotation_line + line + '\n' 
            else:
                output_msa += line + '\n'
        else:
            output_msa += '\n'
    ff = open(output_msa_path, 'wb')
    ff.write(output_msa)
    ff.close()
    return None
