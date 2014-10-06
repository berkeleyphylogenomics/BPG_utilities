#!/usr/bin/env python

import os, glob, re, shlex, subprocess

hmm_name_re = re.compile('bpg(\d+)\.hmm')

def append_family_hmm(library_f, book_dir):
  filenames = os.listdir(book_dir)
  max_accession = 0
  hmm_name_with_max_accession = ''
  for filename in filenames:
    m = hmm_name_re.search(filename)
    if m:
      accession = int(m.group(1))
      if accession > max_accession:
        max_accession = accession
        hmm_name_with_max_accession = filename
  if max_accession > 0:
    full_filename = os.path.join(book_dir, hmm_name_with_max_accession)
    n = subprocess.Popen(['cat', full_filename],
                          stdout = library_f)
  
def main():
  dirs_with_newbooks = glob.glob("/clusterfs/ohana/bpg/pfacts/bpg0/bpg01*")
  library_name = '/clusterfs/ohana/bpg/pfacts/phylofacts_family_hmms'
  f = open(library_name, "w")
  for dir in dirs_with_newbooks:
    bookdirs = os.listdir(dir)
    os.chdir(dir)
    for bookdir in bookdirs:
      append_family_hmm(f, bookdir)
  f.close()
  os.chdir('/clusterfs/ohana/bpg/pfacts')
  n = subprocess.call(['hmmpress', 'phylofacts_family_hmms'])

if __name__ == '__main__':
  main()
