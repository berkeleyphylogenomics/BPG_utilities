#!/usr/bin/env python

import os, sys, re, string
from optparse import OptionParser

def main():
  usage = "%prog [options] orthologs_csv_file taxon1 taxon2 [taxon3] ..."
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) < 3:
    opt_parser.error("Incorrect number of arguments")
  csv_file = args[0]
  if not os.path.exists(csv_file):
    opt_parser.error("Couldn't find %s" % csv_file)
  f = open(csv_file)
  header_line = f.readline()
  lines = f.readlines()
  f.close()
  print "%d lines in %s" % (len(lines), csv_file)
  taxon_re = re.compile("|".join([string.upper(arg) for arg in args[1:]]))
  indices_to_print = set()
  indices_to_print.add(0)
  indices_to_print.add(1)
  headers = header_line.split(',')
  for i in xrange(0,len(headers)):
    m = taxon_re.match(headers[i])
    if m != None:
      indices_to_print.add(i)
  print ','.join([headers[i] for i in indices_to_print])
  for line in lines:
    fields = line.rstrip().split(',')
    m = taxon_re.match(fields[1].split('_')[1])
    if m != None:
      print ','.join([fields[i] for i in indices_to_print])

if __name__ == '__main__':
  main()
