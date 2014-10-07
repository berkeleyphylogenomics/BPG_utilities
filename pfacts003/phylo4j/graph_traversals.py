'''
This file contains functions that traverse the graph and return objects about its structure. 

'''
from pfacts003.phylo4j.consts import *

def traverse(node, max_depth, max_nodes, rel_types, c_depth, c_nodes, v_nodes, c_edges):
    # recursive traversal - this should be rewritten as a cypher query, this is
    # really inefficient...
    curr_nodes = dict(c_nodes)
    curr_edges = dict(c_edges)
    v_nodes.append(node._id)
    if len(curr_nodes) > max_nodes:
        # here we will return a message saying this object is too big
        raise Exception
    if c_depth == max_depth:
        return ({},{})
    for rel in node.match():
        if rel.type in rel_types:
            # add these nodes
            if rel.start_node._id not in curr_nodes:
                p = rel.start_node.get_properties()
                if 'nodeType' in p:
                    nt = p['nodeType'].split('.')[-1]
                else:
                    nt = 'N/A'
                curr_nodes[rel.start_node._id] = {'id':rel.start_node._id, 'type': nt, 'color':'blue'} #, 'type':rel.start_node.type}
            if rel.end_node._id not in curr_nodes:
                p = rel.end_node.get_properties()
                if 'nodeType' in p:
                    nt = p['nodeType'].split('.')[-1]
                else:
                    nt = 'N/A'
                curr_nodes[rel.end_node._id] = {'id':rel.end_node._id, 'type': nt,
                                                'color':'blue'} #, 'type':rel.end_node.type}
            if node == rel.start_node:
                next_node = rel.end_node # outgoing edge
            else:
                next_node = rel.start_node # incoming edge
            # put this edge in the list
            if rel.start_node._id not in curr_edges:
                curr_edges[rel.start_node._id] = {rel.end_node._id : {'type':rel.type, 'id':rel._id}}
            elif rel.end_node._id not in curr_edges[rel.start_node._id]:
                curr_edges[rel.start_node._id][rel.end_node._id] = {'type':rel.type, 'id':rel._id}
            if next_node._id not in v_nodes:
                (n_nodes, n_edges) = traverse(next_node, max_depth, max_nodes, rel_types, c_depth + 1, curr_nodes.items(), v_nodes, curr_edges.items())
                curr_nodes = dict(curr_nodes.items() + n_nodes.items())
                curr_edges = dict(curr_edges.items() + n_edges.items())
    return (curr_nodes, curr_edges)

def explorer_topology_data(node, request_get):
    try:
        # try to get the parameters
        getstring_dict = dict(request_get.lists())
        relationship_types = getstring_dict['relationship_types[]']
        depth = int(getstring_dict['depth'][0])
    except:
        return {}

    try:
        (nodes, edges) = traverse(node, depth, EXPLORER_MAX_NODES, relationship_types, 0, [], [], [])
    except:
        return {'error': 'This subgraph contains too many nodes to display, please adjust the allowed edge types to the right.'}
    return {'nodes': nodes, 'edges': edges} 
