'''
This is a module of commands used in structure determinatio.
Author: Curt Hansen
Created: June 7, 2012
Modified:
'''


import os,sys,math
import pf_connect as db
from Bio.PDB import *
import math_functions as m


def convert_aa_3char_to_1char(aa):
    aadict = {'GLU':'E', 'LYS':'K', 'ARG':'R', 'PHE':'F', 'ASP':'D', 'LEU':'L', 'SER':'S', 'ILE':'I', 'THR':'T',
              'VAL':'V', 'GLY':'G', 'HIS':'H', 'ALA':'A', 'ASN':'N', 'GLN':'Q', 'PRO':'P', 'MET':'M', 'TRP':'W',
              'TYR':'Y', 'CYS':'C'}
    return aadict[aa]


def get_list_residues(chain,form):
    if form not in ['i','n']:
        sys.exit("Form variable must equal 'i' or 'n'!")
    results = []
    for residue in chain:
        if is_aa(residue):
            if form == 'i': #indices
                results.append(residue.id[1])
            else:
                results.append(convert_aa_3char_to_1char(residue.get_resname()))
    return results


def comp_calpha_distance(res1,res2):
    if not (is_aa(res1) and is_aa(res2)):
        sys.exit("Both inputs must be BioPython PDB residues!")
    return m.comp_l2_norm(res1['CA'].get_coord()-res2['CA'].get_coord())


def get_list_of_chains(dbCur,protID):
    results = db.get_query_list(dbCur,"SELECT p.chain_id, u.from_residue, u.to_residue FROM uniprot_pdb_chain as u, \
      pdb_chain as p, uniprot as un WHERE u.pdb_chain_id=p.id AND un.id=u.uniprot_id AND un.accession = %s \
      ORDER BY p.chain_id;",[protID])
    return results


def print_PDB_file_contents(pdbfile):
    parser = PDBParser()
    structure = parser.get_structure('X',pdbfile)
    print "Analyzing contents of",pdbfile
    print "Structure has %s model(s)" % len(structure)
    for model in structure:
        print "Model %s has %s chain(s)" % (model.id, len(model))
        for chain in model:
            print "Chain %s has %s residue(s)" % (chain.id, len(chain))
            for residue in chain:
                print "Residue %s has %s atom(s)" % (residue.get_resname(), len(residue))
                for atom in residue:
#                    print "Atom %s has coordinates %s" % (atom.get_name(), atom.get_vector())
                    print "Atom %s " % (atom.get_name())
