#!/usr/bin/env python2.7

import os
import sys
import time
from rdflib.graph import plugin, Graph, URIRef, Namespace
from rdflib.sparql import SPARQLGraph

def read_owl_file(filename):
    '''Returns the rdf object from eading an owl file.'''
    g = Graph()
    g.parse(filename)
    count = 0
    bp = Namespace("http://www.biopax.org/release/biopax-level2.owl#")
    relations = list(g.subject_objects(bp["protein"]))
    for x,y in relations:
        count +=1
        if count > 20:
            break
        print x,y
#     for subj, obj, predic in g:
#         count += 1
#         if count > 11:
#             break
#         print subj, obj, predic



def get_sequences(owl_obj):
    '''Getting sequence characters from the owl object.'''
    pass


def main(filename):
    '''Main glue function.'''
    owl_obj = read_owl_file(filename)
    sequences = get_sequences(owl_obj)



if __name__ == "__main__":
    start_time = time.time()
    main('/clusterfs/ohana/external/BioCyc/16.0/biocyc-flatfiles-all-w-ecoli/ecoli/16.0/data/biopax-level2.owl')
    print 'Elapsed time: %s' % (time.time() - start_time)
