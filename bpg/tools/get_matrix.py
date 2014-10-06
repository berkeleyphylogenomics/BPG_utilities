#!/usr/bin/python

import os
import sys
import re
import glob

if len(sys.argv) != 3:
    print 'Usage: program taxon_list_file results_directory'
    sys.exit(1)

# read taxon list #

mammal_list = []

handle = open(sys.argv[1], 'r')

for line in handle:
    line = line.strip()
    if line:
        mammal_list.append(line.replace(' ', '_'))

handle.close()

# print matrix @

print '\t',

for x in mammal_list:
    print '%s\t' % x,

print ''

files = glob.glob(sys.argv[2] + '/*.a2m')

for file in files:
    file = file.strip()

    phog = file.split('/')[len(file.split('/'))-1]
    phog = phog[0:len(phog)-4]

    gene_presence = []

    for i in range(len(mammal_list)):
        gene_presence.append(0)

    handle = open(file, 'r')

    for line in handle:
        if line[0] == '>':
            taxon = line.split('_seqh')[0]
            taxon = taxon[1:len(taxon)]

            for i in range(len(mammal_list)):
                if mammal_list[i] == taxon:
                    gene_presence[i] += 1

    handle.close()

    print '%s\t' % phog,

    for x in gene_presence:
        print '%d\t' % x,

    print ''

