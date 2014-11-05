#!/usr/bin/env python

import os, sys
from optparse import OptionParser

def main():
  if len(sys.argv) < 3:
    print "Usage: %s [-s seed] [-p pdb_id or -c comparative model]" % sys.argv[0]
    sys.exit(0)

  Parser = OptionParser('usage: run_discern_on_seed.py <-s seed> <-p pdb_id> <-c comparative model> <-o offset>')
  Parser.add_option('-s', '--seed_ID',
    dest = 'seed_ID',
    default = True,
    help = "seed ID")
  Parser.add_option('-p', '--pdb_ID',
    dest = 'pdb_ID',
    default = False,
    help = "pdb ID (with chain ID)")
  Parser.add_option('-c', '--model_ID',
    dest = 'model_ID',
    default = False,
    help = "model ID (with chain ID)")
  Parser.add_option('-o', '--offset',
    dest = 'offset',
    default = False,
    help = "structure position offset")
  (CmdLineOps, Args) = Parser.parse_args()
 
  seed = CmdLineOps.seed_ID
  if CmdLineOps.pdb_ID:
    pdb_id = CmdLineOps.pdb_ID
    pdb = pdb_id[0:4]
    chain_id = pdb_id[4:5]
  if CmdLineOps.model_ID:
    model_id =CmdLineOps.model_ID
    pdb = model_id[0:4]
    chain_id = model_id[4:5]
    #os.system ('format_modeller.py %s %s.pdb' % (model_id, pdb))
  if CmdLineOps.offset:
    offs = int(CmdLineOps.offset)
  else:
    offs = 0 

  #intrepid_dir = os.path.join(cur_dir(), seed, 'INTREPID')
  intrepid_dir = os.path.join(os.getcwd(), seed, 'INTREPID')
  print os.getcwd()
  print intrepid_dir
  if not os.path.exists(os.path.join(intrepid_dir, 'output.aux')):
    os.system('run_intrepid_on_seed.py %s' % seed)
  os.chdir(intrepid_dir)
  os.environ['LD_LIBRARY_PATH'] \
    = ':'.join([
        "/clusterfs/ohana/software/BALL-1.2/lib/Linux-x86_64-g++_4.1.2/",
        "/clusterfs/ohana/software/BALL-1.2/contrib/lib",
      ])
  os.environ['PATH'] = \
      '/clusterfs/ohana/software/discern/bin:' + os.environ['PATH']
  os.environ['BALL'] = '/clusterfs/ohana/software/BALL-1.3.1'
  naccessdata = "/clusterfs/ohana/software/naccess2.1.1"
  parameter_dir = '/clusterfs/ohana/software/discern/bin/parameters.train-all'
  os.system ("echo $LD_LIBRARY_PATH")
  if CmdLineOps.pdb_ID:
    os.system("cp /clusterfs/ohana/external/pdb/%s/pdb%s.ent.gz %s.pdb.gz" 
            % (pdb[1:3], pdb, pdb))
    os.system("gunzip %s.pdb.gz" % pdb)
  os.system("mv %s.pdb %s.original.pdb" %(pdb, pdb))

  if not chain_id:
    os.system("/clusterfs/ohana/software/intrepid/scripts/repairPDB %s.original.pdb -chain A -offset %d > %s.pdb" % (pdb, offs, pdb))
    chain_id = 'A'
  else:
    os.system("/clusterfs/ohana/software/intrepid/scripts/repairPDB %s.original.pdb -offset %d > %s.pdb" % (pdb, offs, pdb))
       
  os.system("ln -s output.aux intrepid.aux")
  os.system("dsspcmbi %s.pdb %s.dssp" % (pdb, pdb))
  os.system("naccess -r %s/vdw.radii -s %s/standard.data %s.pdb" 
              % (naccessdata, naccessdata, pdb))
  os.system("run_ligsite.pl %s:%s > ligsite.log" % (pdb, chain_id))
  os.system("pdb2graph.pl %s:%s > %s.graph1" % (pdb, chain_id, pdb))
  os.system("analyze-graph %s.graph1 &> analyze.log" % pdb)
  os.system("mv centrality.txt %s.central" % pdb)
  os.system("integrate.pl %s:%s > %s.features" % (pdb, chain_id, pdb))
  os.system("enlargefeatures.pl %s:%s > %s.final-features" 
              % (pdb, chain_id, pdb))
  os.system("rank.pl %s.final-features " % pdb + \
            "%s/weights.txt " % parameter_dir + \
            "%s/mean.txt " % parameter_dir + \
            "%s/stddev.txt " % parameter_dir + \
            " > %s.discern " % pdb)
  #os.system("rank2discern.pl %s.discern intrepid.aux %s > %s.rank"
  #          % (pdb, chain_id, pdb))

if __name__ == '__main__':
  main()
