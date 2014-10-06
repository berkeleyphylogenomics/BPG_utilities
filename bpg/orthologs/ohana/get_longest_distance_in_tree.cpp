// get_longest_distance_in_tree.cpp
// Author: Ruchira S. Datta
// Copyright (c) 2008, Regents of the University of California
// All rights reserved.
//
// Redistiribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// o Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// o Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// o Neither the name of the University of California, Berkeley nor the names
// of its contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR 
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "get_longest_distance_in_tree.h"
#include "find_orthologs_in_tree.h"

#ifdef VERBOSE
#include <iostream>
#endif

double get_longest_distance_in_tree(bpp::Tree &tree,
        const vector<int> &breadth_first_visit_order) {
  const vector<int> nodeIds = tree.getNodesId();
  vector<int>::const_iterator node_iter;
  map<int, ObjectT> object_of_leaf;
  map<int, AttributeT> attribute_of_node;
  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    if (tree.isLeaf(*node_iter)) {
      // Add 1 to avoid the nullObjectValue 0, which is a valid nodeId
      object_of_leaf[*node_iter] = *node_iter + 1;
    }
    attribute_of_node[*node_iter] = 1;
  }
  map<AttributeT, ObjectT > unique_object_with_attribute;
  set<AttributeT> attributes_with_multiple_objects;
  map<int, map<AttributeT, TreeDistanceInfo<ObjectT,
                                      nullObjectValue> *> *>
    distance_from_leaf_to_node;
  find_orthologs_in_tree(tree, breadth_first_visit_order, object_of_leaf,
                        attribute_of_node, unique_object_with_attribute,
                        attributes_with_multiple_objects,
                        distance_from_leaf_to_node,
                        false, true);
  double longest_distance_through_node, longest_distance = 0.0;
  double next_nearest_distance;
  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    longest_distance_through_node =
    (*(distance_from_leaf_to_node[*node_iter]))
      [1]->getDistanceToNearestObjectWithAttributeValue();
    next_nearest_distance 
      = (*(distance_from_leaf_to_node[*node_iter]))
      [1]->getDistanceToNextNearestObjectWithAttributeValue();
    if (next_nearest_distance != -1.0) {
      longest_distance_through_node += next_nearest_distance;
    }
    if (longest_distance_through_node > longest_distance) {
      longest_distance = longest_distance_through_node;
    }
#ifdef VERBOSE
    cout << "Node: " << *node_iter;
    cout << " Furthest leaf: " << 
    (*(distance_from_leaf_to_node[*node_iter]))
      [1]->getNearestObjectWithAttributeValue();
    cout << " Next furthest leaf: " <<
    (*(distance_from_leaf_to_node[*node_iter]))
      [1]->getNextNearestObjectWithAttributeValue();
    cout << " Span: " << longest_distance_through_node;
    cout << endl;
#endif
  }
  return longest_distance;
}
