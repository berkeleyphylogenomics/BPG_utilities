#!/usr/bin/env python2.6

from optparse import OptionParser
from Bio import AlignIO
import os, cStringIO
from subprocess import  check_call

def main():
  parser = OptionParser(usage='%prog [options] alignment_path')
  parser.add_option('--name', dest='name', default='',
                    help="Name of the HMM to be created")
  (options, args) = parser.parse_args()
  if len(args) < 1:
    parser.error('Must supply the path to an alignment file.')

  alignment_path = args[0]
  basename = options.name
  if basename == '':
    basename = os.path.splitext(os.path.split(alignment_path)[1])[0]

  f = open(alignment_path)
  alignment = AlignIO.read(f, 'fasta')
  f.close()
  # Convert alignment to Stockholm format
  # Save in a string for later use, as we will insert a line
  s = cStringIO.StringIO()
  AlignIO.write([alignment], s, "stockholm")
  stockholm_lines = s.getvalue().split('\n')
  s.close()

  # Get the first record from the alignment
  for record in alignment:
    break

  # Write the Stockholm format file
  stockholm_fname = basename + '.stockholm'
  f = open(stockholm_fname, "w")
  # Write everything before the // line
  f.write('\n'.join(stockholm_lines[:-2]))
  f.write('\n')
  # Write the #=RF line indicating which columns are aligned
  f.write("#=GC RF ")
  for ch in record.seq:
    if ch.isupper() or ch == '-':
      f.write("x")
    else:
      f.write(".")
  f.write("\n")
  # Write the // line
  f.write('\n'.join(stockholm_lines[-2:]))
  f.write('\n')
  f.close()

  # create initial HMMER3 hmm

  hmmfname = basename + '.hmm'
  cmd = 'hmmbuild --cpu 8 --amino  --hand %s %s  &> hmmbuild.out' \
        % (hmmfname, stockholm_fname)
  os.system(cmd)

if __name__ == '__main__':
  main()
