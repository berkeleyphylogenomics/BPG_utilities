#!/usr/bin/python

import sys
import os
import re
import MySQLdb
from Bio.Nexus.Trees import Tree

# usage
if len(sys.argv) < 3:
    print 'usage: %s M species_filename' % sys.argv[0]
    print '  where M is the minimum number of matching taxa'
    print '          for an alignment to be returned, and'
    print '        species_filename is the name of a file'
    print '          with one species name per line.'
    print '%s finds ortholog groups containing at least' % sys.argv[0]
    print 'M of the given species names and returns the alignments and'
    print 'trees for those taxa.'
    sys.exit(1)

# get min species to include
M = int(sys.argv[1])

# get user species names
spnames = []

handle = open(sys.argv[2])

for line in handle:
    line = line.strip()
    
    if len(line) > 0:
        spnames.append(line)

handle.close()

# Connect to MySQL.

dbname = 'pfacts002'
dbhost = 'phylodb'
user   = 'pfacts002'
passwd = 'fritzLang'

db  = MySQLdb.connect(db=dbname, host=dbhost, user=user, passwd=passwd)
cur = db.cursor(MySQLdb.cursors.DictCursor)

# get taxon id for species name (common or scientific)
sql = 'SELECT id FROM ncbi_taxonomy WHERE scientific_name in (\'%s\'' % spnames[0]

for i in range(1,len(spnames)):
    sql += ', \'%s\'' % spnames[i]

sql += ') or common_name in (\'%s\'' % spnames[0]

for i in range(1,len(spnames)):
    sql += ', \'%s\'' % spnames[0]

sql += ')'

n_rows = cur.execute(sql)

if n_rows == 0:
    print 'BIG OL\' ERROR FOR... YOU! I didn\'t find any valid species names in your list'
    sys.exit(0)

taxonids = []

for i in range(n_rows):
    taxonids.append(cur.fetchone()['id'])

sql  = 'SELECT tree_node.id, tree_node.tree_id, superorthologous_node_id'
sql += ' FROM tree_node, protein_sequence, genbank_uniprot'
sql += ' WHERE protein_sequence_id <> 0'
sql += ' AND protein_sequence_id = protein_sequence.id'
sql += ' AND genbank_uniprot_id = genbank_uniprot.id'
sql += ' AND ncbi_taxid in (%s' % taxonids[0]

for i in range(1,len(taxonids)):
    sql += ', %s' % taxonids[i]

sql += ')'

n_rows = cur.execute(sql)

if n_rows == 0:
    print 'OOOH, OUCH, no orthologs found'
    sys.exit(0)

# collapse into superortholog nodes
superortho_nodes = set([])

for i in range(n_rows):
    row = cur.fetchone()

    supernode = int(row['superorthologous_node_id'])

    if supernode != 0:
        superortho_nodes.add((row['tree_id'], row['superorthologous_node_id']))

## warning, real work below ##

phog_bpg_dict = {}
alignments    = {}

# get taxa for each superortho-node
for (tree, supernode) in superortho_nodes:
    sql  = 'SELECT id, left_id, right_id FROM tree_node WHERE id=%s' % supernode
    cur.execute(sql)
    row = cur.fetchone()

    myltid = row['left_id']
    myrtid = row['right_id']

    sql  = 'SELECT book.scopid_or_bpgid, ncbi_taxonomy.id, ncbi_taxonomy.scientific_name, sequence.id'
    sql += ' FROM tree, book, alignment, tree_node, protein_sequence, sequence, genbank_uniprot, ncbi_taxonomy'
    sql += ' WHERE tree.id = %s' % tree
    sql += ' AND tree_node.left_id >= %s' % myltid
    sql += ' AND tree_node.right_id <= %s' % myrtid
    sql += ' AND protein_sequence_id <> 0'
    sql += ' AND protein_sequence_id = protein_sequence.id'
    sql += ' AND protein_sequence.genbank_uniprot_id = genbank_uniprot.id'
    sql += ' AND tree_node.tree_id = tree.id'
    sql += ' AND ncbi_taxid = ncbi_taxonomy.id'
    sql += ' AND tree.alignment_id = alignment.id'
    sql += ' AND alignment.book_id = book.id'
    sql += ' AND sequence_id = sequence.id'

    n_rows = cur.execute(sql)

    alignmentid = ''
    phogid = supernode;
    species_inparalog_list = {}

    for i in range(n_rows):
        row = cur.fetchone()

        alignmentid = row['scopid_or_bpgid']

        if row['id'] in taxonids:
            myname = row['scientific_name']

            if myname not in species_inparalog_list:
                species_inparalog_list[myname] = set([])
                
            species_inparalog_list[myname].add(row['sequence.id'])

    # only try to add if at least M species are present
    ntaxa = len(species_inparalog_list)
    
    if ntaxa >= M:
        phog_bpg_dict[phogid] = alignmentid        

        # check to make sure none of the same genes are present
        good_to_add = True
        todelete = set([])

        for (phogacc, aln) in alignments.iteritems():
            for (species, genelst) in aln.iteritems():
                for gene in genelst:
                    for (k,v) in species_inparalog_list.iteritems():
                        if gene in v:
                            existing_len = len(aln)

                            if existing_len >= ntaxa:
                                good_to_add = False
                                break

                            else:
                                todelete.add(phogacc)
    
        if good_to_add:
            alignments[phogid] = species_inparalog_list

        for x in todelete:
            del alignments[x]




### print alignment data ###
for (phog, aln) in alignments.iteritems():
    print '%s (%s) [%d]' % (phog, phog_bpg_dict[phog], len(aln))

    for (sp, genelist) in aln.iteritems():
        print '%s [' % sp.replace(' ', '_'),
        for gene in genelist:
            print gene,

        print ']'
    
    print ''





# done getting info from the DB, now time to look up the alignments, etc

resultdir = '%s.%d.results' % (sys.argv[2], M)

if not os.path.exists(resultdir):
    os.mkdir(resultdir)

for phogacc in alignments:
    aln = alignments[phogacc]
    bpg = phog_bpg_dict[phogacc]

    oldid_newid = {}

    for (sp, lst) in aln.iteritems():
        for gene in lst:
            oldid_newid[gene] = '%s_bpgseq%s' % (sp.replace(' ', '_'), gene)

    if bpg[0:3] == 'bpg':
        dirbase = '/home/bpg/pfacts/'

        ### get the tree ###

        trefname = dirbase +  bpg[0:6] + '/' + bpg[0:9] + '/user/' + bpg[0:9] + '.nj'
        handle = open(trefname, 'r')

        treestr = ''

        for line in handle:
            treestr += line.strip()

        handle.close()

        for oldid in oldid_newid:
            treestr = treestr.replace('bpgseq%d' % oldid, oldid_newid[oldid], 1)

        mytreeobj = Tree(tree=treestr, rooted=True)

        ### prune taxa we don't want ###

        alltaxa = mytreeobj.get_taxa()
        badtaxa = []
        slowest_inparalogs = {}

        for taxon in alltaxa:
            if taxon not in oldid_newid.values():
                badtaxa.append(taxon)

            else:
                sp = taxon.split('_bpgseq')[0]

                if sp in slowest_inparalogs:
                    (old_taxon, old_brlen) = slowest_inparalogs[sp]
                    new_brlen = mytreeobj.node(mytreeobj.search_taxon(taxon)).get_data().branchlength

                    if new_brlen < old_brlen:
                        slowest_inparalogs[sp] = (taxon, mytreeobj.node(mytreeobj.search_taxon(taxon)).get_data().branchlength)
                        badtaxa.append(old_taxon)

                    else:
                        badtaxa.append(taxon)

                else:
                    slowest_inparalogs[sp] = (taxon, mytreeobj.node(mytreeobj.search_taxon(taxon)).get_data().branchlength)

        for taxon in badtaxa:
            mytreeobj.prune(taxon)
        
        ### print tree file ###

        treestr = mytreeobj.to_string(branchlengths_only=True, plain_newick=True)
        
        outfname = '%s/%d.tre' % (resultdir, phogacc)
        outf = open(outfname, 'w')
        print >>outf, treestr
        outf.close()



        ### now get the alignment ###

        alnfname = dirbase + bpg[0:6] + '/' + bpg[0:9] + '/user/shmms/' + bpg[0:9] + '.subfam'
        outfname = '%s/%d.a2m' % (resultdir, phogacc)
        outfile  = open(outfname, 'w')
        handle   = open(alnfname, 'r')

        doneit = set([])
        line = handle.readline()

        while line:
            if line[0] == '>':
                line = line.strip()
                line = line[7:len(line)]
                seqid = int(line)

                if seqid in oldid_newid and seqid not in doneit:
                    doneit.add(seqid)
                    saveit = False

                    for k in slowest_inparalogs:
                        (theid, brlen) = slowest_inparalogs[k]
                        
                        if theid == oldid_newid[seqid]:
                            saveit = True
                            break

                    if saveit:
                        myseq = ''

                        line = handle.readline()

                        while line and line[0] != '>':
                            line   = line.strip()
                            line   = line.replace('.', '')
                            line   = re.sub(r'[a-z]+', '', line)
                            myseq += line
                            line   = handle.readline()

                        print >>outfile, '>%s' % oldid_newid[seqid]
                        print >>outfile, myseq

                    else:
                        line = handle.readline()

                else:
                    line = handle.readline()

            else:
                line = handle.readline()

        handle.close()
        outfile.close()

"""
        alnfname = '/home/bpg/pfacts/'

        alnfname += bpgacc[0:6] + '/' + bpgacc[0:9] + '/user/shmms/' + bpgacc[0:9] + '.subfam'

        outfname = '%s/%d.a2m' % (resultdir, alignments_phogacc[bpgacc])
        outf = open(outfname, 'w')

        handle = open(alnfname, 'r')

        line = handle.readline()

        doneit = set([])

        while line:
            if line[0] == '>':
                myseqid = int(line[7:len(line)].strip())

                found_it = False

                if myseqid not in doneit:

                    for (k,v) in alignments[bpgacc].iteritems():
                        if v == myseqid:
                            found_it = True
                            myheader = '%s_bpgseq%s' % (k.replace(' ', '_'), v)

                            # read seq
                            line = handle.readline()

                            mysequence = ''

                            while line and line[0] != '>' and line[0] != '%':
                                nextl = line.strip().replace('.', '')
                                nextl = re.sub(r'[a-z]+', '', nextl) 
                                mysequence += nextl
                                line = handle.readline()

                            print >>outf, '>%s' % myheader
                            print >>outf, mysequence

                            doneit.add(myseqid)

                            oldid_newid[myseqid] = myheader

                            break

                if not found_it:
                    line = handle.readline()

            else:
                line = handle.readline()

        handle.close()
        outf.close()

        # remove cols with all gaps
        cmd = 'removeGappyColumns 99 %s > xyz_temp.txt.fa' % outfname
        os.system(cmd)
        os.system('mv xyz_temp.txt.fa %s' % outfname)

        # now get the (sub)tree

        trefname = '/home/bpg/pfacts/' +  bpgacc[0:6] + '/' + bpgacc[0:9] + '/user/' + bpgacc[0:9] + '.nj'
        handle = open(trefname, 'r')

        treestr = ''

        for line in handle:
            treestr += line.strip()

        handle.close()

        for oldid in oldid_newid:
            treestr = treestr.replace('bpgseq%d' % oldid, oldid_newid[oldid], 1)

        mytreeobj = Tree(tree=treestr, rooted=True)

        alltaxa = mytreeobj.get_taxa()
        badtaxa = []

        for taxon in alltaxa:
            if taxon not in oldid_newid.values():
                badtaxa.append(taxon)

        for taxon in badtaxa:
            mytreeobj.prune(taxon)

        treestr = mytreeobj.to_string(branchlengths_only=True, plain_newick=True)

        outfname = '%s/%d.tre' % (resultdir, alignments_phogacc[bpgacc])
        outf = open(outfname, 'w')
        print >>outf, treestr
        outf.close()

    else:
        pass
        #print 'no bpg: %s' % bpgacc

for (k,v) in alignments.iteritems():
    print 'ALIGNMENT: %s [phog %s]' % (k, alignments_phogacc[k])

    for (x,y) in v.iteritems():
        print '  %s %s' % (x.replace(' ', '_'), y)

    print ''

"""
