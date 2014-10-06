#!/usr/bin/python
import sys
from pfacts003.phylofacts.models import *

if  len(sys.argv) < 2:
    print "usage: specify family accession"
    sys.exit(1)

accession=sys.argv[1][4:]
family=Family.objects.get(id=accession)
canonical_tree_node=family.canonical_root_node()
c_tree_method=str(canonical_tree_node.tree.method)
seqs=family.canonical_root_node().get_included_leaves()
trees=Tree.objects.filter(family=family)
methods=[str(tree.method) for tree in trees]
methods.sort()
has_ml_tree=has_nj_tree='NULL'
if methods[0]=='ml'and methods[1] =='nj':
    has_ml_tree='y'
    has_nj_tree='y'
elif methods[0]=='ml':
    has_ml_tree='y'
elif methods[0]=='nj':
    has_nj_tree=='y';
print "bpg0%s,%s,%s,%s,%s" %(accession,c_tree_method,has_ml_tree,has_nj_tree,len(seqs))
