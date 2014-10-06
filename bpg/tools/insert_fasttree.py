#!/usr/bin/env python

import os, re
from optparse import OptionParser
from bpg.makebook.insertFamilyIntoDB import insertMLTree

def main():
  parser = OptionParser(usage='%prog [options] alignment_path family_accession')
  (options, args) = parser.parse_args()
  if len(args) != 2:
    parser.error('Must supply alignment_path and family_accession')
  if not os.path.exists(args[0]):
    parser.error('Alignment_path %s not found' % args[0])
  alignment_path = args[0]
  accession_re = re.compile('bpg\d\d\d\d\d\d\d')
  m = accession_re.match(args[1])
  if not m or len(m.group(0)) != len(args[1]):
    parser.error('%s is not a valid BPG family accession' % args[1])
  family_accession = args[1]
  if insertMLTree(alignment_path, family_accession):
    print "Success: ML tree was inserted."
  else:
    print "Error: ML tree was NOT inserted!"

if __name__ == '__main__':
  main()
