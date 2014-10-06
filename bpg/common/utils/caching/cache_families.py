#!/usr/bin/python

import subprocess as sp
import sys
from sys import argv

from pfacts003.phylofacts.models import *
num_families = 0

if len(argv) < 3:
    print '''
    Error: you need to provide the starting family id (omit the
    starting bpg0 part) and the number of families to cache.  Try
    again'''
    sys.exit(1)
try:
    id = int(argv[1])
    num_families_to_cache = int(argv[2])
except:
    print "Error: both the family id and the number of families should be 'intable'!"
    sys.exit(1)

families = Family.objects.order_by("id")

print "num_families = %d" % num_families 
print "num_families_to_cache = %d" % num_families_to_cache
print "starting family = %s" % families[0]
print "ending family = %s" % families[int(num_families_to_cache) - 1]

for family in families[id : id + int(num_families_to_cache)]:
    treenode = family.canonical_root_node()
    if family.status != "bad" and family.active:
        print family
        try:
            fc = FamilyCache.objects.get(family_id = family.id)
            fc.cached_data = treenode.data_for_caching()
            fc.save()
        except:
            fc = FamilyCache(family_id = family.id, 
                    cached_data = treenode.data_for_caching())
            fc.save()
