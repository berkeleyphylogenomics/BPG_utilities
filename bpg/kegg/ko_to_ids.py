#!/usr/bin/python

import os, sys, string

def write_genes(file_handle, organism, gene_line):
  for gene in gene_line[17:].rstrip().split():
    if gene.find('(') >= 0:
      gene = gene.split('(')[0]
    file_handle.write("%s:%s\n" % (organism, gene))

f = open("/usr/local/blastdb/kegg/genes/ko")
lines = f.readlines()
f.close()
os.chdir("/usr/local/blastdb/kegg/genes/ko_ids")

inGenes = False
print "%d lines in ko" % len(lines)
for line in lines:
  if line[0:5] == 'ENTRY':
    entry = line.split()[1]
    inGenes = False
  elif line[0:5] == 'GENES':
    inGenes = True
    id_handle = open("%s.id" % entry, "w")
    organism = string.lower(line[12:15])
    write_genes(id_handle, organism, line)
  elif inGenes: 
    if line[0] == ' ':
      if line[12] != ' ':
        organism = string.lower(line[12:15])
      write_genes(id_handle, organism, line)
    else:
      id_handle.close()
      inGenes = False
