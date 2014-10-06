#! /usr/bin/python

import os
import sys
import re
from Bio import PDB
import psycopg2
import psycopg2.extras

# 
# Input: PDB ID p:q (p=pdb id, q=chain id).
# # This directory must contain the input files: p.pdb

def run(runstr, cmndname):
    cmd = runstr
    
    if debug > 0:
        print cmd

    else:
        print 'running %s...' % cmndname,
        sys.stdout.flush()

    returned = os.system(cmd)

    if debug > 0:
        print 'returned=%d' % returned
        sys.stdout.flush()

    elif returned == 0:
        print 'done.'
        sys.stdout.flush()

    else:
        print 'oops, %s FAILED.' % cmndname
        sys.stdout.flush()


bindir = "/clusterfs/ohana/software/discern/bin"
naccessdata = "/clusterfs/ohana/software/naccess2.1.1"
debug = 0

if len(sys.argv) < 5:
    print 'program <pdb id> <weight file> <mean parameter file> <stddev parameter file>'
    sys.exit(1)

pdbid      = sys.argv[1]
weightfile = sys.argv[2]
meanfile   = sys.argv[3]
stddevfile = sys.argv[4]
bias       = -6.339

j = pdbid.split(':')[0]
k = pdbid.split(':')[1]


# extract sequence from pdb file to seqfile.fa
p = PDB.PDBParser()
structure = p.get_structure('X', '%s.pdb' % j)
chain     = structure[0][k.upper()]

ppb = PDB.PPBuilder()
seqs = ppb.build_peptides(structure)

handle = open('seqfile.fa', 'w')
print >>handle, '>%s:%s' % (j, k)
print >>handle, seqs[0].get_sequence()
handle.close()


run ("%s/dsspcmbi %s.pdb %s.dssp" % (bindir, j, j), 'dssp')
run ("%s/naccess -r %s/vdw.radii  -s %s/standard.data  %s.pdb" % (bindir, naccessdata, naccessdata, j), 'naccess')
run ("%s/run_ligsite.pl > ligsite.log %s" % (bindir, pdbid), 'ligsite')

# run intrepid
run ("intrepid.py seqfile.fa > intrepid.log; cp intrepid.output.aux intrepid.aux", 'intrepid')

# Generates a graph file from a PDB file
run ("%s/pdb2graph.pl %s:%s > %s.graph1" % (bindir, j, k, j), 'pdb2graph')

# Generates centrality features in centrality.txt
run ("%s/analyze-graph %s.graph1 &> analyze.log; mv centrality.txt %s.central" % (bindir, j, j), 'analyze-graph')

# Generates a per-residue features file 
run ("%s/integrate.pl %s:%s > %s.features" % (bindir, j, k, j), 'integrate')

# Combines features from neighboring residues
# Also outpus neighbors.txt and distance.txt that contains lists of neighbors and their distances
run ("%s/enlargefeatures.pl %s:%s > %s.final-features" % (bindir, j, k, j), 'enlarge features')


run ("%s/rank.pl %s.final-features %s %s %s %s > %s.rank" % (bindir, j, weightfile, meanfile, stddevfile, bias, j), 'rank')

# put results (in j.rank) into the db, with method='discern' #

# get family identifier from intrepid output #

'''handle = open('intrepid-alignment.familyid', 'r')
familyid = handle.readline().strip()
handle.close()

### db connection constants ###

dbname = 'pfacts003_test'
user = 'webuser'
password = 'w3zx19Ko'

### end db connection constants  ###

# connect to db

conn = psycopg2.connect("dbname='%s' user='%s' host='db' password='%s'" % (dbname, user, password))
cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

handle = open('%s.rank' % j, 'r')

patt = re.compile('(\d+)\s+(.+)')

for line in handle:
    line = line.strip()

    match = patt.match(line)

    residue = match.group(1)
    score   = match.group(2)

    cur.execute('INSERT INTO residue_score (family_id, residue_number, score, method) VALUES (%s, %s, %s, \'discern\')' % (familyid, residue, score))

handle.close()

conn.commit()

print 'family: %s' % familyid
'''
print 'discern finished.'

