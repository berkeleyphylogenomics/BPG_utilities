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

typedef pair<string, vector<string> *>
  SequenceHeadersForTaxon;
typedef map<int, SequenceHeadersForTaxon>
  SequenceHeadersForSingleTaxonMap;

void find_inparalogs_in_tree(bpp::Tree * tree) {
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
  SequenceHeadersForSingleTaxonMap
    inparalogs_under_node;
  vector<string> *sequence_headers;
  vector<int> nodeIds = tree->getNodesId();
  set<string>::iterator taxon_iter;
  string empty = "";
  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    sequence_headers = new vector<string>;
    inparalogs_under_node[*node_iter] 
      = make_pair("", sequence_headers);
  }
  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    const string &sequence_header = sequence_header_of_sequence[*node_iter];
    const string &taxon = taxon_of_sequence_header[sequence_header];
    inparalogs_under_node[*node_iter].first = taxon;
    inparalogs_under_node[
                      *node_iter].second->push_back(sequence_header);
  }
  vector<int> nodes_to_visit;
  int rootId = tree->getRootId();
  nodes_to_visit.push_back(rootId);
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
  int nodeId, parentId;
  vector<string>::iterator inparalog_iter;
  for (int i = nodes_to_visit.size() - 1; i >= 1; --i) {
    nodeId = nodes_to_visit[i];
    parentId = tree->getFatherId(nodeId);
    if (inparalogs_under_node[nodeId].first != "") {
      if (inparalogs_under_node[parentId].first == "" ||
          inparalogs_under_node[parentId].first ==
          inparalogs_under_node[nodeId].first) {
        inparalogs_under_node[parentId].first
          = inparalogs_under_node[nodeId].first;
        for (inparalog_iter = inparalogs_under_node[nodeId].second->begin();
            inparalog_iter != inparalogs_under_node[nodeId].second->end();
            ++inparalog_iter) {
          inparalogs_under_node[parentId].second->push_back(*inparalog_iter);
        }
      } else {
        inparalogs_under_node[parentId].first = "";
        inparalogs_under_node[parentId].second->clear();
      }
    } else {
      inparalogs_under_node[parentId].first = "";
      inparalogs_under_node[parentId].second->clear();
    }
  }
  if (inparalogs_under_node[rootId].first != "" &&
      inparalogs_under_node[rootId].second->size() > 1) {
    for (inparalog_iter = inparalogs_under_node[rootId].second->begin();
          inparalog_iter != inparalogs_under_node[rootId].second->end();
          ++inparalog_iter) {
      cout << *inparalog_iter << ",";
    }
    cout << endl;
  } else {
    for (int i = 1; i < nodes_to_visit.size(); i++) {
      nodeId = nodes_to_visit[i];
      parentId = tree->getFatherId(nodeId);
      if (inparalogs_under_node[nodeId].first != "" &&
          inparalogs_under_node[parentId].first == "" &&
          inparalogs_under_node[nodeId].second->size() > 1) {
        for (inparalog_iter = inparalogs_under_node[nodeId].second->begin();
              inparalog_iter != inparalogs_under_node[nodeId].second->end();
              ++inparalog_iter) {
          cout << *inparalog_iter << ",";
        }
        cout << endl;
      }
    }
  }
  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    delete inparalogs_under_node[*node_iter].second;
  }
  return;
}
