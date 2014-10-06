#include "find_orthologs_in_uniprot_tree.h"
#include <set>
#include <string>
#include <iostream>
#include <boost/config.hpp>
#include <boost/regex.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

typedef map<string, pair<double, string> > 
  DistanceAndSequenceHeaderForTaxon;
typedef map<int, DistanceAndSequenceHeaderForTaxon *> 
  TaxonToNodeDistanceAndSequenceHeaderMap;

void find_orthologs_in_tree(bpp::Tree * tree) {
  vector<int> leafIds = tree->getLeavesId();
  map<string, int> tree_node_id_of_sequence_id;
  regex uniprot_patt("[^\\|]+\\|[^|]+\\|([^\\|]+)");
  smatch seq_matches, taxon_matches;
  vector<int>::iterator node_iter;
  regex taxon_patt("[^_]+_([^_]+)");
  map<int, string> sequence_header_of_sequence;
  map<string, string> taxon_of_sequence_header;
  set<string> taxa;
  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    if (regex_match(tree->getNodeName(*node_iter), seq_matches, uniprot_patt)){
      string uniprot_id = "";
      uniprot_id += seq_matches[1];
      tree_node_id_of_sequence_id[uniprot_id] = *node_iter;
      sequence_header_of_sequence[*node_iter] = uniprot_id;
      if (regex_match(uniprot_id, taxon_matches, taxon_patt)){
        taxon_of_sequence_header[uniprot_id] = taxon_matches[1];
        taxa.insert(taxon_matches[1]);
      }
    }
  }
  TaxonToNodeDistanceAndSequenceHeaderMap 
    distance_and_sequence_header_from_taxon_to_node;
  vector<int> nodeIds = tree->getNodesId();
  set<string>::iterator taxon_iter;
  string empty = "";
  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    distance_and_sequence_header_from_taxon_to_node[*node_iter]
      = new DistanceAndSequenceHeaderForTaxon;
    for (taxon_iter = taxa.begin(); taxon_iter != taxa.end(); ++taxon_iter) {
     (*(distance_and_sequence_header_from_taxon_to_node[*node_iter]))
        [*taxon_iter] = make_pair(-1.0, empty);
    }
  }
  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    const string &sequence_header = sequence_header_of_sequence[*node_iter];
    const string &taxon = taxon_of_sequence_header[sequence_header];
    (*(distance_and_sequence_header_from_taxon_to_node[*node_iter]))
      [taxon].first = 0.0;
    (*(distance_and_sequence_header_from_taxon_to_node[*node_iter]))
      [taxon].second = sequence_header_of_sequence[*node_iter];
  }
  vector<int> nodes_to_visit;
  nodes_to_visit.push_back(tree->getRootId());
  for(int i = 0; i < nodes_to_visit.size(); i++) {
    vector<int> children = tree->getSonsId(nodes_to_visit[i]);
    for (node_iter = children.begin(); node_iter != children.end();
          ++node_iter) {
      nodes_to_visit.push_back(*node_iter);
    }
  }
  if (nodes_to_visit.size() != tree->getNumberOfNodes()) {
    cout << "Error: did not push all nodes onto visit stack" << endl;
  }
  double taxon_to_node_distance, taxon_to_parent_distance_through_node,
          current_taxon_to_parent_distance;
  int nodeId, parentId;
  for (int i = nodes_to_visit.size() - 1; i >= 1; --i) {
    nodeId = nodes_to_visit[i];
    parentId = tree->getFatherId(nodeId);
    for (taxon_iter = taxa.begin(); taxon_iter != taxa.end(); ++taxon_iter) {
      taxon_to_node_distance
        = (*(distance_and_sequence_header_from_taxon_to_node[nodeId]))
          [*taxon_iter].first;
      taxon_to_parent_distance_through_node 
        = taxon_to_node_distance + tree->getDistanceToFather(nodeId);
      current_taxon_to_parent_distance
        = (*(distance_and_sequence_header_from_taxon_to_node[parentId]))
          [*taxon_iter].first;
      if (taxon_to_node_distance >= 0.0 && 
          (current_taxon_to_parent_distance == -1.0 ||
            taxon_to_parent_distance_through_node 
              < current_taxon_to_parent_distance)) {
        (*(distance_and_sequence_header_from_taxon_to_node[parentId]))
          [*taxon_iter].first = taxon_to_parent_distance_through_node;
        (*(distance_and_sequence_header_from_taxon_to_node[parentId]))
          [*taxon_iter].second =
        (*(distance_and_sequence_header_from_taxon_to_node[nodeId]))
          [*taxon_iter].second;
      }
    }
  }
  int childId;
  double taxon_to_child_distance_through_node; 
  double current_taxon_to_child_distance;
  for (int i = 0; i < nodes_to_visit.size(); i++) {
    nodeId = nodes_to_visit[i];
    vector<int> children = tree->getSonsId(nodes_to_visit[i]);
    for (int j = 0; j < children.size(); j++) {
      childId = children[j];
      for (taxon_iter = taxa.begin(); 
            taxon_iter != taxa.end(); ++taxon_iter) {
        taxon_to_node_distance
          = (*(distance_and_sequence_header_from_taxon_to_node[nodeId]))
            [*taxon_iter].first;
        taxon_to_child_distance_through_node 
          = taxon_to_node_distance + tree->getDistanceToFather(childId);
        current_taxon_to_child_distance
          = (*(distance_and_sequence_header_from_taxon_to_node[childId]))
            [*taxon_iter].first;
        if (taxon_to_node_distance >= 0.0 && 
            (current_taxon_to_child_distance == -1.0 ||
              taxon_to_child_distance_through_node 
                < current_taxon_to_child_distance)) {
          (*(distance_and_sequence_header_from_taxon_to_node[childId]))
            [*taxon_iter].first = taxon_to_child_distance_through_node;
          (*(distance_and_sequence_header_from_taxon_to_node[childId]))
            [*taxon_iter].second =
          (*(distance_and_sequence_header_from_taxon_to_node[nodeId]))
            [*taxon_iter].second;
        }
      }
    }
  }
  typedef adjacency_list <vecS, vecS, undirectedS> Graph;
  Graph G;
  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    const string &sequence_header = sequence_header_of_sequence[*node_iter];
    const string &own_taxon = taxon_of_sequence_header[sequence_header];
    for (taxon_iter = taxa.begin(); taxon_iter != taxa.end(); 
          ++taxon_iter) {
      const string &ortholog =
        (*(distance_and_sequence_header_from_taxon_to_node[*node_iter]))
          [*taxon_iter].second;
      int ortholog_node = tree_node_id_of_sequence_id[ortholog];
      if (ortholog_node > *node_iter &&
            (*(distance_and_sequence_header_from_taxon_to_node[ortholog_node]))
            [own_taxon].second == sequence_header) {
          add_edge(*node_iter, ortholog_node, G);
      }
    }
  }
  vector<int> component(num_vertices(G));
  int num_components = connected_components(G, &component[0]);
  map<int, vector<int> *> nodes_in_components;
  map<int, vector<int> *>::iterator nodes_in_components_iter;
  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    nodes_in_components_iter = nodes_in_components.find(component[*node_iter]);
    if (nodes_in_components_iter == nodes_in_components.end()) {
      nodes_in_components[component[*node_iter]] = new vector<int>;
    }
    nodes_in_components[component[*node_iter]]->push_back(*node_iter);
  }
  cout << "Component,";
  cout << "Sequence,";
  for (taxon_iter = taxa.begin(); taxon_iter != taxa.end(); ++taxon_iter) {
    cout << *taxon_iter << "_ortholog,";
    cout << *taxon_iter << "_ortholog_is_reciprocal,";
  }
  printf("\n");
  for (nodes_in_components_iter = nodes_in_components.begin();
        nodes_in_components_iter != nodes_in_components.end();
        ++nodes_in_components_iter) {
    for (node_iter = nodes_in_components_iter->second->begin(); 
          node_iter != nodes_in_components_iter->second->end(); ++node_iter) {
      const string &sequence_header = sequence_header_of_sequence[*node_iter];
      cout << nodes_in_components_iter->first << ",";
      cout << sequence_header << ",";
      const string &own_taxon = taxon_of_sequence_header[sequence_header];
      for (taxon_iter = taxa.begin(); taxon_iter != taxa.end(); 
            ++taxon_iter) {
        const string &ortholog =
          (*(distance_and_sequence_header_from_taxon_to_node[*node_iter]))
            [*taxon_iter].second;
        cout << ortholog << ",";
        int ortholog_node = tree_node_id_of_sequence_id[ortholog];
        if ((*(distance_and_sequence_header_from_taxon_to_node[ortholog_node]))
            [own_taxon].second == sequence_header) {
          cout << "Y,";
        } else {
          cout << "N,";
        }
      }
      cout << endl;
    }
  }
      
  for (nodes_in_components_iter = nodes_in_components.begin();
        nodes_in_components_iter != nodes_in_components.end();
        ++nodes_in_components_iter) {
    delete nodes_in_components_iter->second;
  }
  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    delete distance_and_sequence_header_from_taxon_to_node[*node_iter];
  }
  return;
}
