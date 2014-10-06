#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser
from pfacts003.phylofacts.models import SequenceHeader
import cStringIO

seqhdr_re = re.compile('(SEQHDR[0-9]*)')
capturing_seqhdr_re = re.compile('(SEQHDR([0-9]*))')

def replace_sequence_header_numbers_with_identifiers_in_tree(tree_string):
  output = cStringIO.StringIO()
  tree_tokens = seqhdr_re.split(tree_string)
  for i in range(len(tree_tokens)):
    if i % 2 == 0:
      output.write(tree_tokens[i])
    else:
      m = capturing_seqhdr_re.match(tree_tokens[i])
      sequence_header = SequenceHeader.objects.get(id = m.group(2))
      if sequence_header.uniprot:
        output.write(sequence_header.uniprot.uniprot_identifier)
      else:
        output.write(sequence_header.identifier())
  return output.getvalue()

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-m", "--method", dest="method", default="nj",
                      help="tree method for which to find PHOG-T's.")
  opt_parser.add_option("-v", "--verbose", dest="verbose",
              action="store_true", default=True,
              help="Whether to print verbose output")
  opt_parser.add_option("-q", "--quiet", dest="verbose",
              action="store_false", default=True,
              help="Whether to suppress verbose output")
  (options, args) = opt_parser.parse_args()
  if len(args) < 1:
    opt_parser.error("Incorrect number of arguments")
  else:
    bpg_accession = args[0]
  verbose = options.verbose
  tree_path = os.path.join('/clusterfs/ohana/bpg/pfacts/',
                            bpg_accession[0:4],
                            bpg_accession[0:7],
                            bpg_accession,
                            "%s.%s" % (bpg_accession, options.method))
  f = open(tree_path)
  tree_string = f.read()
  f.close()
  translated_string \
    = replace_sequence_header_numbers_with_identifiers_in_tree(tree_string)
  sys.stdout.write(translated_string)
  sys.stdout.write('\n')

if __name__ == '__main__':
  main()
