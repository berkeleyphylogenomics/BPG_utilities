#!/usr/bin/env python

import os, sys, subprocess
from Bio import SeqIO

def main():
  indir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5'
  os.chdir(indir)
  f = open("non_uniprot_accessions.err")
  lines = f.readlines()
  f.close()

  ensembl_ids = set()
  uniprot_accessions_of_ensembl_id = {}
  for line in lines:
    ensembl_id = line.strip().split(':')[2].split()[0]
    ensembl_ids.add(ensembl_id)
    uniprot_accessions_of_ensembl_id[ensembl_id] = set()

  inf = open("9606_homo_sapiens.fasta")
  human_records = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))
  inf.close()

  inf = open("10090_mus_musculus.fasta")
  mouse_records = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))
  inf.close()

  inf = open("7955_danio_rerio.fasta")
  zebrafish_records = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))
  inf.close()

  outf = open("sequences_without_uniprot_accessions.fa", "w")
  records = set()
  for ensembl_id in ensembl_ids:
    key = 'ENSEMBL:%s' % ensembl_id
    if ensembl_id[0:4] == 'ENSP':
      records.add(human_records[key])
    elif ensembl_id[0:7] == 'ENSMUSP':
      records.add(mouse_records[key])
    elif ensembl_id[0:7] == 'ENSDARP':
      records.add(zebrafish_records[key])
  SeqIO.write(list(records), outf, "fasta")
  outf.close()

  p = subprocess.Popen(['blastall',
                '-p', 'blastp',
                '-i', 'sequences_without_uniprot_accessions.fa',
                '-o', 'sequences_without_uniprot_accessions.blast_results',
                '-m' '9',
                '-e', '1e-10',
                '-b', '20',
                '-v', '20',
                '-F', 'F',
                '-d', 'UniProt/current/protein',
                ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, error = p.communicate(input=None)

if __name__ == '__main__':
  main()
