#!/usr/bin/env python

import os, sys, glob
from optparse import OptionParser
from Bio import SeqIO

Parser = OptionParser("%prog [options] <pfam_accession>")
Parser.add_option('-o', dest = 'file', default = '', 
  help = "send removed subfamilies and bpg accessions to FILE")
Parser.add_option('-q', '--quiet', action = 'store_false', dest = 'verbose',
  default = True, help = "do not print progress messages")
Parser.add_option('--test', action = 'store_true', dest = 'test',
  default = False, help = "run unit testing suite and exit")
(CmdLineOps, Args) = Parser.parse_args()
verbose = CmdLineOps.verbose

if CmdLineOps.test:
  import unittest
  #Run the tests, but there are none for now.
  sys.exit(0)

starting_directory = os.getcwd()

def main():
  if len(Args) != 1:
    Parser.print_usage()
    sys.exit(2)

  pfam_accession = Args[0]
  pfam_dirs = glob.glob('/clusterfs/ohana/bpg/pfam/%s.*' % pfam_accession)
  if len(pfam_dirs) == 0:
    if verbose: print "No PhyloFacts families for %s" % pfam_accession
    sys.exit(0)
  os.chdir(pfam_dirs[0])
  subfam_dirs = [dir for dir in glob.glob('subfam*') if os.path.isdir(dir)]

  if len(subfam_dirs) < 2:
    if verbose: print "We don't have multiple PhyloFacts families for %s, so nothing to do"\
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
      if verbose: print "Subfam %s does not appear to have been inserted, skipping" \
        % subfam
    elif len(families) > 1:
      if verbose: print "Subfam %s appears to contain bad families, skipping" % subfam
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
  accessions_already_covered = set()
  newly_covered_by_max = accessions_in_subfam[max_covering_subfam]
  if verbose: print "Covering subfams"
  while len(accessions_to_cover) > 0:
    print "%s, %s, %g, %g, %g" % (max_covering_subfam, 
                    bpg_accession_of_subfam[max_covering_subfam],
                    num_genes_covered_by_max_covering_subfam,
                    len(accessions_in_subfam[max_covering_subfam] \
                    & accessions_already_covered),
                    len(accessions_in_subfam[max_covering_subfam]))
#    print "  ", ','.join(list(newly_covered_by_max))
    accessions_already_covered = accessions_already_covered.union\
        (accessions_in_subfam[max_covering_subfam])
    for accession in accessions_in_subfam[max_covering_subfam]:
      line = [accession, str(accession in accessions_to_cover)]
      previous_subfams = []
      if accession not in accessions_to_cover:
        for subfam in subfams_used:
          if accession in accessions_in_subfam[subfam]:
            previous_subfams.append(subfam)
      previous_subfams = '\t'.join(previous_subfams)
      line.append(previous_subfams)
      print '\t'.join(line)

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
    if CmdLineOps.file:
      os.chdir(starting_directory)
      if verbose: print "Writing redundant subfamilies to file %s"\
        % CmdLineOps.file
      redundancy_file = open(CmdLineOps.file, 'a')
      for subfam in accessions_in_subfam.keys():
        if subfam not in subfams_used:
          redundancy_file.write("%s,%s\n"\
            % (subfam, bpg_accession_of_subfam[subfam]))
      redundancy_file.close()

if __name__ == '__main__':
  main()
