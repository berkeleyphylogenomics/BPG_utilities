#!/usr/bin/python

import os
import sys
sys.path = ['/home/ruchira/ohana_repository'] + sys.path
from bpg.common.utils.dir_of_family import get_dir_of_family_accession
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] family_accession_file"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-i", "--input_file", dest="input", 
                        help="family_accession_file")
  opt_parser.add_option("-d", "--description", dest="description",
                        default='090811',
                        help="The description for download files")

  (options, args) = opt_parser.parse_args()
  if not options.input:
    opt_parser.error('Please specify the family_accession_file')
  if not os.path.exists(options.input):
    opt_parser.error("Couldn't find %s" % options.input)
  else:
    input_file = options.input
  f = open(input_file, "rU")
  family_accessions = [line.rstrip() for line in f.readlines()]
  f.close()
  hmm_file = 'downloads/%s.hmms' \
              % (options.description)
  print "Writing %s with HMMs" % (hmm_file)
  for family_accession in family_accessions:
      dir = get_dir_of_family_accession(family_accession)
      hmm = os.path.join(dir, '%s.hmm' % family_accession)
      try:

# put newly added hmms to a new hmm file (for download)
        cmd ="cat %s >> %s" % (hmm,hmm_file)
        os.system (cmd)
      except IOError:
        print "Could not write %s to the hmm file" % hmm

if __name__ == '__main__':
  main()
