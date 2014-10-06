// find_inparalogs_in_tree.cpp
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


#include "find_inparalogs_in_tree.h"
AttributeT nullAttributeValue = 0;

//template<class AttributeT, AttributeT nullAttributeValue>
void find_inparalogs_in_tree(const bpp::Tree &tree, 
                              const vector<int> &breadth_first_visit_order,
                              map<int, AttributeT > &attribute_of_node) {
  // DO NOT clear attribute_of_node here, as it already has the attributes of
  // the leaves.
  typedef map<int, AttributeT > AttributeNodeMap;
//  typename 
  AttributeNodeMap::const_iterator attribute_of_node_iter;
  int nodeId, childId;
  bool isInparalogous;
  for(int i = breadth_first_visit_order.size() - 1; i >= 0; --i) {
    nodeId = breadth_first_visit_order[i];
    if (tree.isLeaf(nodeId)) {
      continue;
    }
    const vector<int> &children = tree.getSonsId(nodeId);
    isInparalogous = false;
    if (children.size() > 0) {
      childId = children[0];
      attribute_of_node_iter = attribute_of_node.find(childId);
      if (attribute_of_node_iter != attribute_of_node.end()) {
        attribute_of_node[nodeId] = attribute_of_node[childId];
        isInparalogous = true;
        for(int j = 1; j < children.size() && isInparalogous; ++j) {
          childId = children[j];
          attribute_of_node_iter = attribute_of_node.find(childId);
          if (attribute_of_node_iter == attribute_of_node.end() ||
              attribute_of_node[nodeId] != nullAttributeValue &&
              attribute_of_node[childId] != nullAttributeValue &&
              attribute_of_node[childId] != attribute_of_node[nodeId]) {
            attribute_of_node.erase(nodeId);
            isInparalogous = false;
          } else if (attribute_of_node[nodeId] == nullAttributeValue) {
            attribute_of_node[nodeId] = attribute_of_node[childId];
          }
        }
      }
    }
  }
}
