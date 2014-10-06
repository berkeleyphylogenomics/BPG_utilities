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

  unaligned_seqs_file = '%s.seqs.fa' % bpg_accession
  blastable_database_index = '%s.pin' % unaligned_seqs_file

  if not os.path.exists(blastable_database_index):
    print 'No BLASTable database found.  Will create.'
    starttime = time.time()
    print 'Creating file of unaligned sequences from %s.' % bpg_accession
    sys.stdout.flush()
    p = subprocess.Popen(['alignedfasta_to_unalignedfasta.py', 
                        '%s.a2m' % bpg_accession,
                        '%s.seqs.fa' % bpg_accession], stderr=subprocess.PIPE)
    status = os.waitpid(p.pid, 0)[1]
    
    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)
    starttime = time.time()
    print 'Creating BLASTable database of sequences from %s.' % bpg_accession
    sys.stdout.flush()
    p = subprocess.Popen(['formatdb', '-o', 'T', '-l', 'formatdb.log', '-i',
                        '%s.seqs.fa' % bpg_accession], stderr=subprocess.PIPE)
    status = os.waitpid(p.pid, 0)[1]
    
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

  tree_objects = Tree.objects.filter(family__id = family_id, method='ml')
  if not tree_objects:
    starttime = time.time()
    print 'Inserting ML tree...',
    sys.stdout.flush()

    status = insertMLTree('final_trimmed_mafft_ungapped.afa', bpg_accession)

    if not status:
      sys.stderr.write("Error: ML tree was NOT inserted!\n")

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)

  if not os.path.exists('ml_tree_ordered.a2m'):
    print 'Tree-ordered alignment file not found, creating.'
    family = Family.objects.get(id = family_id)

    aligned_seq_of_seq_header = {}
    tree_node_alignment_objects = TreeNodeAlignment.objects.filter(
                  tree_node = family.canonical_root_node())
    starttime = time.time()
    print 'Reading alignment from database...',
    sys.stdout.flush()
    for object in tree_node_alignment_objects:
      aligned_seq_of_seq_header[object.sequence_header] \
          = object.aligned_sequence
    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)
    tree = Tree.objects.get(family__id = family_id, method='ml')
    starttime = time.time()
    print 'Reading order from database...',
    sys.stdout.flush()
    leaves = TreeNode.objects.filter(tree = tree, sequence_header__isnull =
                                  False).order_by('left_id')
    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)
    starttime = time.time()
    print 'Writing reordered alignment...',
    sys.stdout.flush()
    f = open('ml_tree_ordered.a2m','w')
    for leaf in leaves:
      f.write('>%s\n' % leaf.sequence_header.header)
      f.write('%s\n' % aligned_seq_of_seq_header[leaf.sequence_header].chars)
    f.close()
    cmd = 'rm %s.a2m; ln -s ml_tree_ordered.a2m %s.a2m' %  (bpg_accession,
                                                            bpg_accession)
    p = subprocess.Popen(cmd, shell=True)
    status = os.waitpid(p.pid, 0)[1]
    starttime = time.time()
    print 'done...',
    sys.stdout.flush()

  tree_file = bpg_accession + '.ml'

  if os.path.exists(tree_file):
    starttime = time.time()
    print 'Running find_orthologs on ML tree...',
    sys.stdout.flush()

    outf = open('%s.find_orthologs2.ml.out' % bpg_accession, 'w')
    errf = open('%s.find_orthologs2.ml.err' % bpg_accession, 'w')
    cmd = 'find_orthologs --book %s --method ml' % bpg_accession
    print cmd

    p = subprocess.Popen(shlex.split(cmd), stdout=outf, stderr=errf)

    status = os.waitpid(p.pid, 0)[1]

    errf.close()
    outf.close()

    endtime = time.time()
    print 'done. %s' % getTimeStr(starttime, endtime)
  
    p = subprocess.Popen(['touch', 'run_phog2_ml_again_done'] )
    status = os.waitpid(p.pid, 0)[1]
  else:
    print "Didn't find %s, skipped." % tree_file

  overall_endtime = time.time()
  print "Done polishing family %s. Overall: %s" \
      % (bpg_accession, getTimeStr(overall_starttime, overall_endtime))
  
  p = subprocess.Popen(['touch', 'polish_family_20110430_done'] )
  status = os.waitpid(p.pid, 0)[1]

if __name__ == '__main__':
  main()
