#!/usr/bin/env python

import os, subprocess, shlex, time, glob, sys
from optparse import OptionParser
from bpg.makebook.buildFamily import getTimeStr
from bpg.makebook.insertFamilyIntoDB import insertMLTree

def main():
  # parse command line options
  usage = "%prog [options] bpg_accession"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if len(args[0]) != 10 or args[0][0:3] != 'bpg':
    opt_parser.error('Argument must be a bpg accession like bpg0123456')
  bpg_accession = args[0]
  try:
    family_id = int(bpg_accession[3:])
  except ValueError:
    opt_parser.error('Argument must be a bpg accession like bpg0123456')
  family_dir = '/clusterfs/ohana/bpg/pfacts/%s/%s/%s' % (bpg_accession[0:4],
                                                          bpg_accession[0:7],
                                                          bpg_accession) 
  if not os.path.exists(family_dir):
    opt_parser.error('Family %s not found on the filesystem.' % bpg_accession)

  print "Polishing family %s in %s" % (bpg_accession, family_dir)
  os.chdir(family_dir)

  overall_starttime = time.time()

  if not os.path.exists(bpg_accession + '.alignmentconservation.csv'):
    starttime = time.time()
    print 'computing alignment conservation and inserting...',
    sys.stdout.flush() 

    errf = open('compute_alignment_conservation.err','w')
    p = subprocess.Popen(['compute_alignment_conservation.py', bpg_accession],
                          stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

  f = open(bpg_accession + '.hmm')
  line = f.readline().strip()
  f.close()

  if line != 'HMMER3/b [3.0 | March 2010]':
    starttime = time.time()
    print "Old HMMER3 version line: ", line
    print 'Rebuilding the HMM using the production version of HMMER3...',
    sys.stdout.flush()

    stockholm_fname = glob.glob('*.stockholm')[0]
    outf = open('hmmer3build.out', 'w')
    errf = open('hmmer3build.err', 'w')
    cmd = 'hmm3build -n %s --hand %s.hmm %s' \
          % (bpg_accession, bpg_accession, stockholm_fname)
    p = subprocess.Popen(shlex.split(cmd), stdout=outf, stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    outf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

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

    starttime = time.time()
    print 'Inserting ML tree...',
    sys.stdout.flush()

    status = insertMLTree('final_trimmed_mafft_ungapped.afa', bpg_accession)

    if not status:
      sys.stderr.write("Error: ML tree was NOT inserted!\n")

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    starttime = time.time()
    print 'Running find_orthologs...',
    sys.stdout.flush()

    outf = open('%s.find_orthologs.ml.out' % bpg_accession, 'w')
    errf = open('%s.find_orthologs.ml.err' % bpg_accession, 'w')
    cmd = 'find_orthologs --book %s --method ml' % bpg_accession

    p = subprocess.Popen(shlex.split(cmd), stdout=outf, stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    outf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

    starttime = time.time()
    print 'Running find_thresholded_phogs_django.py...',
    sys.stdout.flush()

    outf = open('%s.find_thresholded_phogs.ml.out' % bpg_accession, 'w')
    errf = open('%s.find_thresholded_phogs.ml.err' % bpg_accession, 'w')
    cmd = 'find_thresholded_phogs_django.py -m ml %s' % bpg_accession

    p = subprocess.Popen(shlex.split(cmd), stdout=outf, stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    outf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

  
  overall_endtime = time.time()
  print "Done polishing family %s. Overall: %s" \
      % (bpg_accession, getTimeStr(overall_starttime, overall_endtime))
  
  p = subprocess.Popen(['touch', 'polish_family_done'] )
  status = os.waitpid(p.pid, 0)[1]

if __name__ == '__main__':
  main()
