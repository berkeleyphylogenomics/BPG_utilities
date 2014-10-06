// write_phogs_at_threshold.cpp
// Author: Ruchira S. Datta
// Copyright (c) 2011, Regents of the University of California
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

#include "write_phogs_at_threshold.h"
#include "find_duplication_distances.h"
#include <iostream>
#include <assert.h>

void printLeavesDescendingFromNode(const bpp::Tree &tree, int nodeId) {
  if (tree.isLeaf(nodeId)) {
    cout << "  " << tree.getNodeName(nodeId) << endl;
  } else {
    const vector<int> &children = tree.getSonsId(nodeId);
    for (int j = 0; j < children.size(); ++j) {
      printLeavesDescendingFromNode(tree, children[j]);
    }
  }
}

void writePHOGsAtThreshold(const bpp::Tree &tree,
    const LeftRightIdsOfNodeIdMap &idMap,
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *> &pair_of_nearest_objects_to_maximal_node,
    const map<int,  double> &greatest_distance_of_maximal_descendant,
    double threshold) {
  map<int, double>::const_iterator max_desc_dist_iter;
  map<int, ObjObjInfo *>::const_iterator int_objobj_iter;
  double duplication_distance;
  int leftId;
  // Write PHOGs satisfying
  // greatest_distance_of_maximal_descendant < threshold <= duplication_distance
  for (max_desc_dist_iter = greatest_distance_of_maximal_descendant.begin();
      max_desc_dist_iter != greatest_distance_of_maximal_descendant.end();
      ++max_desc_dist_iter) {
    if (max_desc_dist_iter->second >= threshold) {
      continue;
    }
    leftId = idMap.find(max_desc_dist_iter->first)->second.first;
    int_objobj_iter 
      = pair_of_nearest_objects_to_maximal_node.find(max_desc_dist_iter->first);
    if (int_objobj_iter == pair_of_nearest_objects_to_maximal_node.end()) {
      // This had better be the root
      assert(leftId == 1);
      duplication_distance = 5000.0;
    } else {
      duplication_distance = getMaximalNodeDistance(max_desc_dist_iter->first,
                            duplication_node_of_objects,
                            pair_of_objects_yielding_duplication_node_distance,
                            pair_of_nearest_objects_to_maximal_node);

    }
    if (threshold <= duplication_distance) {
      cout << "PHOG " << leftId << " ";
      if (max_desc_dist_iter->second <= 0.0) {
        cout << "[0.0,";
      } else {
        cout << "(" << max_desc_dist_iter->second << ",";
      }
      cout << duplication_distance << "]" << endl;
      printLeavesDescendingFromNode(tree, max_desc_dist_iter->first);
    }
  }
}
