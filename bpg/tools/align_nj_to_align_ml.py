#!/usr/bin/env python

import os, glob, sys

def main():
  aln_nj_files = glob.glob('*.aln.nj')
  if len(aln_nj_files) == 0:
    print "No .aln.nj file found, exiting..."
    sys.exit(0)
  else:
    for aln_nj_file in aln_nj_files:
      basename = os.path.splitext(os.path.splitext(aln_nj_file)[0])[0]
      aln_ml_file = '%s.aln.ml' % basename
      if os.path.exists(aln_ml_file):
        print "%s already exists, skipping..." % aln_ml_file
        continue
      f = open(aln_nj_file)
      lines = f.readlines()
      f.close()
      if len(lines) == 0:
        print "%s is empty, skipping..." % aln_nj_file
        continue
      # Expect alignment lines like this:
      # SEQ367 --ALRSIIV--
      fields = lines[0].split()
      if len(fields) != 2:
        print "%s is not in the expected format, skipping..." % aln_nj_file
        continue
      f = open(aln_ml_file, "w")
      # Print the number of sequences and the number of 
      f.write('%d %d\n' % (len(lines), len(fields[1])))
      for line in lines:
        f.write(line)
      f.close()

if __name__ == '__main__':
  main()
