#!/usr/bin/python

import os, sys

"""
Input: PDB coordinate file and optional chain ID.
Output: fasta sequence file for chain, or for chain A if no chain specified.
  Fasta header contains pdb ID + chain ID (e.g. 1mk4A)
Sequence derived from ATOM records.

MOTIVATION: Ruchira says that sometimes the initial methionine is incuded
  in the PDB file but not in the PDB blast DB (or vice versa?)
FATAL DEFICIENCY: usually, some residues are not represented by ATOM records.
  These may even be in the middle of a sequence.
POSSIBLE CORRECTIONS:
  -- use SEQRES records instead. deficiency:  what if we have a file
     (perhaps from non-PDB source) that lacks SEQRES records?
  -- get seqlen from some other record and put X's in for missing residues.
     This seems really bad -- we want the missing residues for sequence
     alignment purposes.
"""

if len(sys.argv) < 2:
  print "Usage: %s <pdb_coord_file> [<chain_ID>]" % sys.argv[0]
  sys.exit(0)

# make dictionary to translate from 3-letter to 1-letter
one_letter_of_three_letter = {}
one_letter_of_three_letter['ACE'] = ''
one_letter_of_three_letter['GLY'] = 'G'
one_letter_of_three_letter['ALA'] = 'A'
one_letter_of_three_letter['VAL'] = 'V'
one_letter_of_three_letter['LEU'] = 'L'
one_letter_of_three_letter['ILE'] = 'I'
one_letter_of_three_letter['SER'] = 'S'
one_letter_of_three_letter['CYS'] = 'C'
one_letter_of_three_letter['THR'] = 'T'
one_letter_of_three_letter['MET'] = 'M'
one_letter_of_three_letter['MSE'] = 'M'
one_letter_of_three_letter['PRO'] = 'P'
one_letter_of_three_letter['PHE'] = 'F'
one_letter_of_three_letter['TRP'] = 'W'
one_letter_of_three_letter['HIS'] = 'H'
one_letter_of_three_letter['LYS'] = 'K'
one_letter_of_three_letter['TYR'] = 'Y'
one_letter_of_three_letter['ARG'] = 'R'
one_letter_of_three_letter['ASP'] = 'D'
one_letter_of_three_letter['GLU'] = 'E'
one_letter_of_three_letter['ASN'] = 'N'
one_letter_of_three_letter['GLN'] = 'Q'
seed_id = sys.argv[1]

pdb_path  = sys.argv[1]
if not os.path.exists(pdb_path):
  print >> sys.stderr, "%s: file %s does not exist." % \
     (sys.argv[0], pdb_path)
  sys.exit(1)

chain = 'A'
if len(sys.argv) > 2:
  chain = sys.argv[2]
  if len(chain) > 1:
    print >> sys.stderr, "%s: chain ID %s should be single character." % \
       (sys.argv[0], chain)
    sys.exit(1)
  chain = chain.upper()

pdb = open(pdb_path, "r")
lines = pdb.read().split('\n')
prev_seqno = -32767
found_chain = False

if lines[0].startswith("HEADER"):
  pdb_id = lines[0][62:66].lower()
  sys.stdout.write('>%s' % pdb_id)
  if chain: sys.stdout.write(chain)
  sys.stdout.write('\n')
else:
  print >> sys.stderr, "%s: pdb file %s does not begin with HEADER record." % \
     (sys.argv[0], pdb_path)
  sys.exit(1)

# Look for ATOM records with out chain ID.
# When the sequence number increments, record a new residue.
for line in lines[1:]:
  if (line.startswith('ATOM') or line.startswith('HETATM')) and \
      (not chain or line[21] == chain):
    found_chain = True
    cur_seqno = int(line[22:26].strip())
    if cur_seqno > prev_seqno:
      sys.stdout.write(one_letter_of_three_letter[line[17:20]])
      prev_seqno = cur_seqno
  elif found_chain and chain and line[0:3] == 'TER':
    break

sys.stdout.write('\n')
