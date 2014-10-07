''' This file contains the api functions for intrepid. '''

import re, tempfile, codecs
from piston.handler import BaseHandler
from datetime import datetime
from pfacts003.phylo4j.consts import *
from pfacts003.phylo4j.table_views import explorer_node_property_table, explorer_relationships_table 
from pfacts003.phylo4j.graph_traversals import explorer_topology_data
from py2neo import neo4j
from django.core.validators import validate_email
from django.core.exceptions import ValidationError
from django.contrib.auth.models import User

class Phylo4jExplorerSearch(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request):
        if 'query' in request.GET:
            q = request.GET['query']
        else:
            return {'error':'s','message':'No query found.'}
        node = None
        # connect to the server
        graph_db = neo4j.GraphDatabaseService(NEO4J_SERVER_ADDRESS)
        # check to see if we have that protein node
        # is it a uniprot accession ?
        if graph_db.get_indexed_node('protein_accession_index','protein_accession_index', q[:6]):
            node = graph_db.get_indexed_node('protein_accession_index','protein_accession_index', q[:6])
            graph_db = None
            return {'node_id': int(node._id)}
        #if graph_db:
            
        # If all else fails, tell the user we couldn't find their query
        graph_db = None
        return {'error':'s','message':'%s was not found in phylo4j.' % (q)}

class Phylo4jExplorerData(BaseHandler):
    allowed_methods = ('GET')

    def read(self, request, id, t):
        graph_db = neo4j.GraphDatabaseService(NEO4J_SERVER_ADDRESS)
        n = graph_db.node(int(id))
        if not n:
            graph_db = None
            return {'error':'s', 'message':'Cannot find node %d' % (int(id))}
        graph_db = None
        if t == 'relationships':
            return explorer_relationships_table(n)
        elif t == 'graph_data':
            return explorer_topology_data(n, request.GET)
        elif t == 'properties':
            # return the properties of this node in a tabular format
            return explorer_node_property_table(n)
        return '<p>Data could not be loaded</p>'
