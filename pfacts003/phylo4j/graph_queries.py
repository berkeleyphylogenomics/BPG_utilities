''' This file contains the cypher queries we use for anything in phylo4j. '''
from pfacts003.phylo4j.consts import *
from py2neo import neo4j

class Interaction:
    # This represents an IntActInteraction (two nodes in bio4j connected by a 
    # PROTEIN_PROTEIN_INTERACTION edge.
    # Make sure to pass the start node correctly, because it assumes that intact_id_1 is the start
    def __init__(self, start_node, end_node, interaction, int_type):
        # what kind of interaction is this?
        if int_type == 'PROTEIN_PROTEIN_INTERACTION':
            self.interaction_source = 'IntAct'
        elif int_type == 'dip_ppi':
            self.interaction_source = 'DIP'
        elif int_type == 'AGREEMENT_SUBTREE_PPI':
            self.interaction_source = 'Interolog'
        else:
            self.interaction_source = 'Unknown'
        self.start_node_accession = start_node['accession']
        self.end_node_accession = end_node['accession']
        try:
            self.num_experiments_supporting = interaction['experiments']
        except:
            pass
        try:
            self.start_node_intact_id = interaction['intact_id_1']
            self.end_node_intact_id = interaction['intact_id_2']
        except:
            pass
        self.start_node_uniprot_identifier = start_node['name']
        self.end_node_uniprot_identifier = end_node['name']
        self.start_node_description = start_node['full_name']
        self.end_node_description = end_node['full_name']
        # ebi made these weird... need to figure out how to link to the intact interaction
        self.interaction_html_link = ''

def get_interactions(query):
    # This returns a list of all the intact interactions for a query, cast into the class above
    graph_db = neo4j.GraphDatabaseService(NEO4J_SERVER_ADDRESS)
    if graph_db.get_indexed_node('protein_accession_index','protein_accession_index', query[:6]):
        this_node = graph_db.get_indexed_node('protein_accession_index','protein_accession_index', query[:6])
    else:
        graph_db = None
        return []

    # Get the Intact interactions from bio4j
    # Current relationship types we have in phylo4j
    interaction_types = ['PROTEIN_PROTEIN_INTERACTION','dip_ppi','AGREEMENT_SUBTREE_PPI']
    
    these_interactions = []
    for t in interaction_types:
        these_interactions += [Interaction(this_node.get_properties(), 
            rel.end_node.get_properties(), rel.get_properties(), t) for rel in 
            graph_db.match(start_node = this_node, rel_type=t)]
    
    # don't know if this is necessary, probably not.
    graph_db = None
    return these_interactions

def get_all_relationships():
    # returns ALL of the possible relationship types in the graph
    graph_db = neo4j.GraphDatabaseService(NEO4J_SERVER_ADDRESS)
    # get the relationship types from the database
    s = [str(x) for x in graph_db.relationship_types]
    # sort them alphabetically
    s.sort()
    # maybe not necessary?
    graph_db = None
    return s
