#!/usr/bin/env python

import os, subprocess, shlex, time, glob, sys
from optparse import OptionParser
from bpg.makebook.buildFamily import getTimeStr
from bpg.makebook.insertFamilyIntoDB import insertMLTree
from pfacts003.phylofacts.models import Family, Tree, TreeNode, \
    TreeNodeAlignment

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if len(args[0]) != 10 or args[0][0:3] != 'bpg':
    opt_parser.error('Argument must be a bpg accession like bpg0164198')
  bpg_accession = args[0]
  try:
    family_id = int(bpg_accession[3:])
  except ValueError:
    opt_parser.error('Argument must be a bpg accession like bpg0164198')
  family_dir = '/clusterfs/ohana/bpg/pfacts/%s/%s/%s' % (bpg_accession[0:4],
                                                          bpg_accession[0:7],
                                                          bpg_accession) 
  if not os.path.exists(family_dir):
    opt_parser.error('Family %s not found on the filesystem.' % bpg_accession)

  print "Polishing family %s in %s" % (bpg_accession, family_dir)
  os.chdir(family_dir)

  overall_starttime = time.time()


  fasttrees = glob.glob('*fasttree*')

  if len(fasttrees) == 0:
    print 'No FastTree ML tree found.  Will create.'
    starttime = time.time()
    print 'Reformatting alignment...',
    sys.stdout.flush()
    errf = open('align_nj_to_align_ml.err', 'w')
    p = subprocess.Popen('align_nj_to_align_ml.py', stderr=errf)
  
    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    
    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    starttime = time.time()
    print 'Running FastTree...',
    sys.stdout.flush()

    unrooted_tree_name = 'final_trimmed_mafft_ungapped.fasttree.ml.tre'
    outf = open(unrooted_tree_name, 'w')
    errf = open('final_trimmed_mafft_ungapped.fasttree.log', 'w')
    p = subprocess.Popen(['FastTree', '-gamma',
                          'final_trimmed_mafft_ungapped.aln.ml'],
                          stdout=outf, stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    outf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    starttime = time.time()
    print 'Rerooting tree at midpoint of the longest span...',
    sys.stdout.flush()

    outf = open('final_trimmed_mafft_ungapped.fasttree.ml.rooted.tre', 'w')
    errf = open('midpoint_reroot.err', 'w')
    cmd = 'export PYTHONPATH=/clusterfs/ohana/software/lib/python2.6; ' \
          + 'midpoint_reroot.py %s ' % unrooted_tree_name
    p = subprocess.Popen(cmd, shell=True, stdout=outf, stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    outf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

  overall_endtime = time.time()
  print "Done polishing family %s. Overall: %s" \
      % (bpg_accession, getTimeStr(overall_starttime, overall_endtime))
  
  p = subprocess.Popen(['touch', 'polish_family_20110103_done'] )
  status = os.waitpid(p.pid, 0)[1]

if __name__ == '__main__':
  main()
