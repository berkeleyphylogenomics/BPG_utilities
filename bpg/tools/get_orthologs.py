#!/usr/bin/python

import sys
import os
import re
import psycopg2
import psycopg2.extras
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

### db connection constants ###

dbname = 'pfacts003_test'
user   = 'bryan'
password = 'mypass'

#user = 'webuser'
#password = 'w3zx19Ko'

### end db connection constants  ###

# connect to db

conn = psycopg2.connect("dbname='%s' user='%s' host='db' password='%s'" % (dbname, user, password))
cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)


# get taxon id for species name (common or scientific)
sql = 'SELECT id FROM ncbi_taxonomy WHERE scientific_name in (\'%s\'' % spnames[0]

for i in range(1,len(spnames)):
    sql += ', \'%s\'' % spnames[i]

sql += ') or common_name in (\'%s\'' % spnames[0]

for i in range(1,len(spnames)):
    sql += ', \'%s\'' % spnames[0]

sql += ')'

cur.execute(sql)

curresults = list(cur.fetchall())
n_rows = len(curresults)

if n_rows == 0:
    print 'BIG OL\' ERROR FOR... YOU! I didn\'t find any valid species names in your list'
    sys.exit(0)

taxonids = []

for i in range(n_rows):
    taxonids.append(curresults[i]['id'])

sql  = 'SELECT tree_node.id, tree_node.tree_id, superorthologous_node_id'
sql += ' FROM tree_node, sequence_header'
sql += ' WHERE tree_node.sequence_id <> 0'
sql += ' AND tree_node.sequence_id = sequence_header.id'
sql += ' AND sequence_header.taxon_id in (%s' % taxonids[0]

for i in range(1,len(taxonids)):
    sql += ', %s' % taxonids[i]

sql += ')'

cur.execute(sql)
curresults = list(cur.fetchall())
n_rows = len(curresults)

if n_rows == 0:
    print 'OOOH, OUCH, no orthologs found'
    sys.exit(0)

# collapse into superortholog nodes
superortho_nodes = set([])

for i in range(n_rows):
    row = curresults[i]

    if row['superorthologous_node_id']:
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

    sql  = 'SELECT family.id, ncbi_taxonomy.id, ncbi_taxonomy.scientific_name, sequence_header.id'
    sql += ' FROM tree, family, tree_node, sequence_header, ncbi_taxonomy'
    sql += ' WHERE tree.id = %s' % tree
    sql += ' AND tree_node.left_id >= %s' % myltid
    sql += ' AND tree_node.right_id <= %s' % myrtid
    sql += ' AND tree_node.sequence_id <> 0'
    sql += ' AND tree_node.sequence_id = sequence_header.id'
    sql += ' AND tree_node.tree_id = tree.id'
    sql += ' AND sequence_header.taxon_id = ncbi_taxonomy.id'
    sql += ' AND tree.family_id = family.id'


#    sql  = 'SELECT book.scopid_or_bpgid, ncbi_taxonomy.id, ncbi_taxonomy.scientific_name, sequence.id'
#    sql += ' FROM tree, book, alignment, tree_node, protein_sequence, sequence, genbank_uniprot, ncbi_taxonomy'
#    sql += ' WHERE tree.id = %s' % tree
#    sql += ' AND tree_node.left_id >= %s' % myltid
#    sql += ' AND tree_node.right_id <= %s' % myrtid
#    sql += ' AND protein_sequence_id <> 0'
#    sql += ' AND protein_sequence_id = protein_sequence.id'
#    sql += ' AND protein_sequence.genbank_uniprot_id = genbank_uniprot.id'
#    sql += ' AND tree_node.tree_id = tree.id'
#    sql += ' AND ncbi_taxid = ncbi_taxonomy.id'
#    sql += ' AND tree.alignment_id = alignment.id'
#    sql += ' AND alignment.book_id = book.id'
#    sql += ' AND sequence_id = sequence.id'

    cur.execute(sql)
    curresults = list(cur.fetchall())
    n_rows = len(curresults)

    alignmentid = ''
    phogid = supernode;
    species_inparalog_list = {}

    for i in range(n_rows):
        row = curresults[i]

        alignmentid = row[0]

        if row[1] in taxonids:
            myname = row[2]

            if myname not in species_inparalog_list:
                species_inparalog_list[myname] = set([])
                
            species_inparalog_list[myname].add(row[3])

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
            oldid_newid[gene] = '%s_seqh%s' % (sp.replace(' ', '_'), gene)

    dirbase = '/home/bpg/pfacts/'

    ### get the tree ###

    sql = 'SELECT tree FROM family_tree WHERE family_id=%s and tree_method=\'ml\'' % bpg
    cur.execute(sql)

    rows = cur.fetchall()

    for row in rows:
        treestr = row['tree']

    for oldid in oldid_newid:
        treestr = treestr.replace('seqh%d' % oldid, oldid_newid[oldid], 1)

    mytreeobj = Tree(tree=treestr, rooted=True)

    ### prune taxa we don't want ###

    alltaxa = mytreeobj.get_taxa()
    badtaxa = []
    slowest_inparalogs = {}

    for taxon in alltaxa:
        if taxon not in oldid_newid.values():
            badtaxa.append(taxon)

        else:
            sp = taxon.split('_seqh')[0]

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

    sql = 'SELECT alignment FROM family_alignment WHERE family_id=%s' % bpg
    cur.execute(sql)

    for row in cur.fetchall():
        alnstr = row['alignment']

    # do some parsing #
    alnarr1 = alnstr.split('\n')

    myalndict = {}

    i = 0
    while i < len(alnarr1):
        if alnarr1[i] and  alnarr1[i][0] == '>':
            oldh  = alnarr1[i].strip()
            seqid = int(oldh[1:len(oldh)])

            if seqid in oldid_newid:
                saveit = False

                for k in slowest_inparalogs:
                    (theid, brlen) = slowest_inparalogs[k]

                    if theid == oldid_newid[seqid]:
                        saveit = True
                        break

                if saveit:
                    myseq = ''

                    i += 1

                    while i < len(alnarr1) and alnarr1[i] and alnarr1[i][0] != '>' and alnarr1[i][0] != '%':
                        line   = alnarr1[i].strip()
                        line   = line.replace('.', '')
                        line   = re.sub(r'[a-z]+', '', line)
                        myseq += line
                        i += 1
            
                    myalndict['>' + oldid_newid[seqid]] = myseq
                
                else:
                    i += 1
            else:
                i += 1
        else:
            i += 1

    outfname = '%s/%d.a2m' % (resultdir, phogacc)
    outfile  = open(outfname, 'w')

    for (header, sequence) in myalndict.iteritems():
        print >>outfile, header
        print >>outfile, sequence

    outfile.close()

