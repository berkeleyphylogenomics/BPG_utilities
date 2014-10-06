#!/usr/bin/env python

from optparse import OptionParser
from pfacts003.phylofacts.models import PDB, PDB_Chain, Sequence
from Bio import SeqIO 
from Bio.SeqUtils import CheckSum

def main():
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()

  f = open("/clusterfs/ohana/external/pdb_rcsb_full", "rU")
  for record in SeqIO.parse(f, "fasta"):
    fields = record.description.split()
    if len(fields) < 2 or fields[1] != 'mol:protein':
      continue
    pdb_id, chain_id = fields[0].split('_')
    pdb_objects = PDB.objects.filter(id__exact = pdb_id)
    if pdb_objects:
      pdb = pdb_objects[0]
    else:
      pdb = PDB.objects.create(id = pdb_id)
    seguid = CheckSum.seguid(record.seq)
    sequence_objects = Sequence.objects.filter(seguid__exact = seguid)
    if sequence_objects:
      sequence = sequence_objects[0]
    else:
      sequence = Sequence.objects.create(chars = record.seq.tostring(),
                                          seguid = seguid)
    pdb_chain_objects = PDB_Chain.objects.filter(pdb__exact = pdb,
                                                chain_id__exact = chain_id)
    if pdb_chain_objects:
      pdb_chain = pdb_chain_objects[0]
      # Update the sequence in case it has changed
      if pdb_chain.full_sequence != sequence:
        pdb_chain.full_sequence = sequence
        pdb_chain.all_residues_have_atom_records_f = None
      if len(fields) >= 4:
        pdb_chain.description = ' '.join(fields[3:])
      pdb_chain.save()
    else:
      if len(fields) >= 4:
        description = ' '.join(fields[3:])
        pdb_chain = PDB_Chain.objects.create(pdb = pdb, chain_id = chain_id,
                                            full_sequence = sequence,
                                            description = description)
      else:
        pdb_chain = PDB_Chain.objects.create(pdb = pdb, chain_id = chain_id,
                                            full_sequence = sequence)
    
if __name__ == '__main__':
  main()
