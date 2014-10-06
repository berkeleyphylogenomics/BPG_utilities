#!/usr/bin/env python

import os, sys, glob
from Bio import SeqIO

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <pfam_accession>" % sys.argv[0]
    sys.exit(0)

  pfam_accession = sys.argv[1]
  pfam_dirs = glob.glob('/clusterfs/ohana/bpg/pfam/%s.*' % pfam_accession)
  if len(pfam_dirs) == 0:
    print "No PhyloFacts families for %s" % pfam_accession
    sys.exit(0)
  os.chdir(pfam_dirs[0])
  subfam_dirs = [dir for dir in glob.glob('subfam*') if os.path.isdir(dir)]

  if len(subfam_dirs) < 2:
    print "We don't have multiple PhyloFacts families for %s, so nothing to do"\
      % pfam_accession
    sys.exit(0)
  
  subfams_containing_accession = {}
  accessions_in_subfam = {}
  bpg_accession_of_subfam = {}

  accessions_to_cover = set()
  num_genes_covered_by_max_covering_subfam = 0
  max_covering_subfam = ""
  for subfam in subfam_dirs:
    os.chdir(subfam)
    families = glob.glob('bpg*.a2m')
    if len(families) < 1:
      print "Subfam %s does not appear to have been inserted, skipping" \
        % subfam
    elif len(families) > 1:
      print "Subfam %s appears to contain bad families, skipping" % subfam
    else:
      accessions_in_subfam[subfam] = set()
      bpg_accession = os.path.splitext(families[0])[0]
      bpg_accession_of_subfam[subfam] = bpg_accession
      f = open(families[0])
      for line in f.readlines():
        if len(line) > 4 and (line[0:4] == '>tr|' or line[0:4] == '>sp|'):
          accession = line[1:].strip().split()[0].split('/')[0]
          accessions_to_cover.add(accession)
          if accession not in subfams_containing_accession:
            subfams_containing_accession[accession] = set()
          subfams_containing_accession[accession].add(subfam)
          accessions_in_subfam[subfam].add(accession)
      f.close()
      num_genes_covered_by_subfam = len(accessions_in_subfam[subfam])
      if num_genes_covered_by_subfam > num_genes_covered_by_max_covering_subfam:
        max_covering_subfam = subfam
        num_genes_covered_by_max_covering_subfam = num_genes_covered_by_subfam
    os.chdir('..')

  subfams_used = set()
  newly_covered_by_max = accessions_in_subfam[max_covering_subfam]
  print "Covering subfams"
  while len(accessions_to_cover) > 0:
    print "%s,%s,%d" % (max_covering_subfam, 
                    bpg_accession_of_subfam[max_covering_subfam],
                    num_genes_covered_by_max_covering_subfam)
    print "  ", ','.join(list(newly_covered_by_max))
    accessions_to_cover -= accessions_in_subfam[max_covering_subfam]
    subfams_used.add(max_covering_subfam)
    max_covering_subfam = ""
    num_genes_covered_by_max_covering_subfam = 0
    for subfam in accessions_in_subfam.keys():
      if subfam not in subfams_used:
        newly_covered = accessions_in_subfam[subfam] & accessions_to_cover
        if len(newly_covered) > num_genes_covered_by_max_covering_subfam:
          max_covering_subfam = subfam
          num_genes_covered_by_max_covering_subfam = len(newly_covered)
          newly_covered_by_max = newly_covered

  if len(accessions_to_cover) > 0:
    print "%d uncovered accessions" % len(accessions_to_cover)
    for accession in accessions_to_cover:
      print accession

  if len(subfams_used) < len(accessions_in_subfam.keys()):
    print "Redundant subfams (can remove):"
    for subfam in accessions_in_subfam.keys():
      if subfam not in subfams_used:
        print "%s,%s" % (subfam, bpg_accession_of_subfam[subfam])

if __name__ == '__main__':
  main()
