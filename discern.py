#!/usr/bin/python

import os
import sys

# open input file and read pdbid:chain
handle = open(sys.argv[1], 'r')
pdbstr = handle.readline().strip()
handle.close()

pdbid   = pdbstr[0:4]
chainid = pdbstr[5]

# get the pdb file from wherever it is
pdbfname = '/clusterfs/ohana/external/pdb_entries/%s/pdb%s.ent' % (pdbid[1:3], pdbid)
os.system('cp %s ./%s.pdb' % (pdbfname, pdbid))

# run discern (web version) on this file
discern = '/clusterfs/ohana/software/discern/bin/discern_web.sh'
os.system('%s %s:%s' % (discern, pdbid, chainid))


