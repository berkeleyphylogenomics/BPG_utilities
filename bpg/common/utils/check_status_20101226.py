#!/usr/bin/env python
import os, subprocess, shlex, glob, sys
from optparse import OptionParser
def main():
  # parse command line options
  usage = "%prog [options]"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  f = open('/home/ruchira/families_20101217.txt')
  bpg_accessions = [line.strip() for line in f.readlines()]
  f.close()
  for bpg_accession in bpg_accessions:
    family_id = int(bpg_accession[3:])
    family_dir = '/clusterfs/ohana/bpg/pfacts/%s/%s/%s' % (bpg_accession[0:4],
                                                            bpg_accession[0:7],
                                                            bpg_accession) 
    if not os.path.exists(family_dir):
      opt_parser.error('Family %s not found on the filesystem.' % bpg_accession)
    print "Checking family %s in %s" % (bpg_accession, family_dir)
    os.chdir(family_dir)
    if os.path.exists('final_trimmed_mafft_ungapped.fasttree.ml.rooted.tre'):
      print 'PASS: %s: rooted FastTree exists' % bpg_accession
      if os.path.exists('%s.ml' % bpg_accession):
        print 'PASS: %s: %s.ml exists' % (bpg_accession, bpg_accession)
      else:
        print 'OOPS: %s: %s.ml does not exist' % (bpg_accession, bpg_accession)
      if os.path.exists('ml_tree_ordered.a2m'):
        print 'PASS: %s: ml_tree_ordered.a2m exists' % bpg_accession
        p = subprocess.Popen(['count_seqs', 
                              'final_trimmed_mafft_ungapped.afa'],
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE)
        output, error = p.communicate(input=None)
        num_seqs = int(output)
        p = subprocess.Popen(['count_seqs', 
                              'ml_tree_ordered.a2m'],
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE)
        output, error = p.communicate(input=None)
        if output and len(output) > 0:
          try:
            num_ml_seqs = int(output)
          except ValueError:
            print 'FAIL: %s: count_seqs ml_tree_ordered.a2m returned %s' \
                % (bpg_accession, output)
          num_ml_seqs = int(output)
          if (num_seqs == num_ml_seqs):
            print 'PASS: %s: alignments have same number of seqs' \
                % bpg_accession
          else:
            print 'FAIL: %s: alignments have different numbers of seqs' \
                % bpg_accession
        else:
          print 'FAIL: %s: count_seqs ml_tree_ordered.a2m returned no output'\
            % bpg_accession
      else:
        print 'FAIL: %s: ml_tree_ordered.a2m does not exist' % bpg_accession
    else:
      print 'INCOMPLETE: %s: rooted FastTree does not exist' % bpg_accession


if __name__ == '__main__':
  main()
