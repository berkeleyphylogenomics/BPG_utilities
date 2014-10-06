#ifndef GRAPH_OF_PHYLOGENETIC_TREE_H
#define GRAPH_OF_PHYLOGENETIC_TREE_H

#include <map>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <Phyl/Tree.h>

using namespace std;
using namespace boost;

typedef enum { bidirectedGT, parentToChildGT, childToParentGT } EdgeConvertT;

struct VertexProperties {
  int tree_node_index;
  string vertex_name;
  bool is_taxon;
};

struct BranchProperties {
  float length;
};

typedef adjacency_list<vecS, vecS, directedS, VertexProperties,
                        BranchProperties > Graph;

typedef graph_traits<Graph>::vertex_descriptor vtx_t;
typedef graph_traits<Graph>::edge_descriptor edge_t;

typedef map<int, vtx_t> TreeNodeVertexMap;

extern pair<Graph *, TreeNodeVertexMap *> 
graph_of_phylogenetic_tree(bpp::Tree * tree);

#endif
