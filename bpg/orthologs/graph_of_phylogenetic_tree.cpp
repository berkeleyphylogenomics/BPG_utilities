#include "graph_of_phylogenetic_tree.h"
#include <iostream>

vtx_t add_vertex_from_treenode(bpp::Tree * tree, int id, Graph &graph,
                                TreeNodeVertexMap &map) 
{
  vtx_t vtx;
  vtx = add_vertex(graph);
  graph[vtx].tree_node_index = id;
  graph[vtx].is_taxon = false;
  map[id] = vtx;
  if (tree->hasNodeName(id)) {
    graph[vtx].vertex_name = tree->getNodeName(id);
  }
  vector<int> children = tree->getSonsId(id);
  vector<int>::iterator i;
  double branch_length;
  vtx_t child_vtx;
  edge_t edge_to_child;
  edge_t edge_from_child;
  for(i = children.begin(); i != children.end(); i++) {
    child_vtx = add_vertex_from_treenode(tree, *i, graph, map);
    branch_length = tree->getDistanceToFather(*i);
    edge_to_child = add_edge(vtx, child_vtx, graph).first;
    graph[edge_to_child].length = branch_length;
    edge_from_child = add_edge(child_vtx, vtx, graph).first;
    graph[edge_from_child].length = branch_length;
  }
  return vtx;
}
    
pair<Graph *, TreeNodeVertexMap *> 
graph_of_phylogenetic_tree(bpp::Tree * tree) {
  Graph * graph = new Graph(0);
  TreeNodeVertexMap * map = new TreeNodeVertexMap();
  add_vertex_from_treenode(tree, tree->getRootId(), *graph, *map);
  return make_pair(graph, map);
}
