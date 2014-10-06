// find_duplication_distances.cpp
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

#include "find_duplication_distances.h"
#include <iostream>

//template<class ObjectT, ObjectT nullObjectValue>
void find_distances_to_maximizing_objects(const bpp::Tree &tree,
        const vector<int> &breadth_first_visit_order,
        const map<ObjectT, map<ObjectT, set<int> *> *> 
          &alternative_nearest_objects,
        const map<ObjectT, int> &leaf_of_object,
        map<int, map<ObjectT, double> *> &distances_to_maximizing_objects) {
  set<ObjectT> maximizing_objects;
  map<ObjectT, map<ObjectT, set<int> *> *>::const_iterator
    obj_objs_nodes_iter;
  map<ObjectT, set<int> *>::const_iterator obj_nodes_iter;
  ObjectT obj0, obj1;
  for (obj_objs_nodes_iter 
        = alternative_nearest_objects.begin();
      obj_objs_nodes_iter
        != alternative_nearest_objects.end();
      ++obj_objs_nodes_iter) {
    obj0 = obj_objs_nodes_iter->first;
    maximizing_objects.insert(obj0);
    for (obj_nodes_iter = obj_objs_nodes_iter->second->begin();
        obj_nodes_iter != obj_objs_nodes_iter->second->end();
        ++obj_nodes_iter) {
      obj1 = obj_nodes_iter->first;
      maximizing_objects.insert(obj1);
    }
  }
  set<ObjectT>::const_iterator max_obj_iter;
  int nodeId;
  for (max_obj_iter = maximizing_objects.begin();
      max_obj_iter != maximizing_objects.end(); ++max_obj_iter) {
    nodeId = leaf_of_object.find(*max_obj_iter)->second;
    distances_to_maximizing_objects[nodeId] = new map<ObjectT, double>;
    (*distances_to_maximizing_objects[nodeId])[*max_obj_iter] = 0.0;
  }
  map<int, map<ObjectT, double> *>::iterator node_obj_dist_iter, 
                                          parent_node_obj_dist_iter;
  map<ObjectT, double>::iterator obj_dist_iter, parent_obj_dist_iter;
  double parent_distance_to_object_through_node,
      parent_distance_to_object,
      node_distance_to_object;
  int parentId;
  ObjectT obj;
  for (int i = breadth_first_visit_order.size() - 1; i >= 1; --i) {
    nodeId = breadth_first_visit_order[i];
    parentId = tree.getFatherId(nodeId);
    node_obj_dist_iter = distances_to_maximizing_objects.find(nodeId);
    if (node_obj_dist_iter == distances_to_maximizing_objects.end()) {
      continue;
    }
    parent_node_obj_dist_iter = distances_to_maximizing_objects.find(parentId);
    if (parent_node_obj_dist_iter 
          == distances_to_maximizing_objects.end()) {
      distances_to_maximizing_objects[parentId] = new map<ObjectT, double>;
    }
    for (obj_dist_iter = distances_to_maximizing_objects[nodeId]->begin();
        obj_dist_iter != distances_to_maximizing_objects[nodeId]->end();
        ++obj_dist_iter) {
      obj = obj_dist_iter->first;
      node_distance_to_object = obj_dist_iter->second;
      parent_distance_to_object_through_node
        = node_distance_to_object + tree.getDistanceToFather(nodeId);
      parent_obj_dist_iter 
        = distances_to_maximizing_objects[parentId]->find(obj);
      if (parent_obj_dist_iter 
          == distances_to_maximizing_objects[parentId]->end()) {
        (*distances_to_maximizing_objects[parentId])[obj]
          = parent_distance_to_object_through_node;
      } else {
        parent_distance_to_object
          = (*distances_to_maximizing_objects[parentId])[obj];
        if (parent_distance_to_object_through_node
            < parent_distance_to_object) {
          (*distances_to_maximizing_objects[parentId])[obj]
            = parent_distance_to_object_through_node;
        }
      }
    }
  }
}

int binary_search(int query, vector<int> &sorted_vector) {
  // We do not use the standard STL binary search algorithm because that merely
  // returns true if the query is in the container and false if it is not.  We
  // expect the query may *not* be in the container.  Instead, this function
  // returns the unique integer i, if any, such that
  // sorted_vector[i] <= query < sorted_vector[i+1].
  // Here sorted_vector is sorted in ascending order.
  int low = 0; int high = sorted_vector.size() - 1;
  if (sorted_vector[low] > query) {
    return -1;
  }
  if (sorted_vector[high] <= query) {
    return high;
  }
  // Now it is the case that low <= high and
  // sorted_vector[low] <= query < sorted_vector[high]
  // which will remain invariant.
  // Since it remains invariant, it cannot be the case at the end of the loop
  // that high = low.  So at the end of the loop, high = low + 1.
  int mid;
  while (high - low > 1) {
    mid = (low + high) / 2; // low <= mid <= high
    if (sorted_vector[mid] <= query) {
      low = mid;  // still sorted_vector[low] = query
    } else {  // query < sorted_vector[mid]
      high = mid; // still query < sorted_vector[high]
    }
  }
  return low;
}

int reverse_binary_search(int query, vector<int> &sorted_vector) {
  // This is like binary_search, but works on a vector in reversed order.
  int low = 0; int high = sorted_vector.size() - 1;
  if (sorted_vector[low] < query) {
    return -1;
  }
  if (sorted_vector[high] >= query) {
    return high;
  }
  // Now it is the case that low <= high and
  // sorted_vector[low] >= query > sorted_vector[high]
  // which will remain invariant.
  // Since it remains invariant, it cannot be the case at the end of the loop
  // that high = low.  So at the end of the loop, high = low + 1.
  int mid;
  while (high - low > 1) {
    mid = (low + high) / 2; // low <= mid <= high
    if (sorted_vector[mid] >= query) {
      low = mid;  // still sorted_vector[low] = query
    } else {  // query > sorted_vector[mid]
      high = mid; // still query > sorted_vector[high]
    }
  }
  return low;
}

double getDuplicationNodeDistance(int nodeId,
  const map<ObjectT, map<ObjectT, DupInfo *> *> 
      &duplication_node_of_objects,
  const map<int, ObjObjInfo *>
          &pair_of_objects_yielding_duplication_node_distance) {
  map<int, ObjObjInfo *>::const_iterator obj_obj_int_iter;
  obj_obj_int_iter 
    = pair_of_objects_yielding_duplication_node_distance.find(nodeId);
  return duplication_node_of_objects.find(
              obj_obj_int_iter->second->getLesserObj())->second->find(
              obj_obj_int_iter->second->getGreaterObj())->second->getHalfSpan();
}

double getMaximalNodeDistance(int nodeId,
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *>
      &pair_of_nearest_objects_to_maximal_node) {
  map<int, ObjObjInfo *>::const_iterator obj_obj_int_iter;
  obj_obj_int_iter = pair_of_nearest_objects_to_maximal_node.find(nodeId);
  // Note that the DupInfo object whose getDuplicationNodeId method is being
  // called here has a getHalfSpan() method, but this does *not* necessarily
  // return the same result as we get by calling getDuplicationNodeDistance().
  // getHalfSpan() returns half the distance between a pair of objects which
  // contributes to making nodeId a maximal node.  But, the
  // duplication_node_distance might be half the distance between some other
  // pair of objects which also got duplicated at that duplication node, and
  // have a greater half span.
  return getDuplicationNodeDistance(
    duplication_node_of_objects.find(
      obj_obj_int_iter->second->getLesserObj())->second->find(
      obj_obj_int_iter->second->getGreaterObj()
                                              )->second->getDuplicationNodeId(),
      duplication_node_of_objects,
      pair_of_objects_yielding_duplication_node_distance);
}

void find_duplication_distances(const bpp::Tree &tree,
        const vector<int> &breadth_first_visit_order,
        const map<int, AttributeT> attribute_of_node,
        const map<int, map<AttributeT, TreeDistanceInfo<ObjectT, 
                                                    nullObjectValue> *> *>
          distance_from_object_with_attribute_to_node,
        const LeftRightIdsOfNodeIdMap &idMap,
        const map<int, int> &nodeIdOfLeftIdMap,
        const map<ObjectT, int> &leaf_of_object,
        const map<int, ObjectT> &object_of_leaf,
        const map<ObjectT, map<ObjectT, set<int> *> *> 
          &alternative_nearest_objects,
        map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
        map<int, ObjObjInfo *>
          &pair_of_objects_yielding_duplication_node_distance,
        map<int, ObjObjInfo *>
          &pair_of_nearest_objects_to_maximal_node,
        map<int,  double> 
          &greatest_distance_of_maximal_descendant) {

  if (alternative_nearest_objects.size() == 0) {
    // There is nothing to do here but set the root
    int rootId = tree.getRootId();
    greatest_distance_of_maximal_descendant[rootId] = -1;
    return;
  }
  // A maximal node is one for which, for some attribute, it is not the case
  // that the node, its parent, and all its siblings have the same nearest
  // object in the tree with that attribute.  There must be at least two
  // distinct objects with the same attribute such that they are the nearest
  // object with that attribute to some of the node, its parent, and its
  // siblings.  All such pairs of objects, and the maximal nodes to which they
  // pertain, have previously been recorded in alternative_nearest_objects.

  // For each pair of objects in alternative_nearest_objects, we will find
  // the most recent common ancestor (MRCA) of that pair of objects.  That is
  // the duplication node for the pair of objects.  We also will find the tree
  // distance between the two objects.  If the tree obeyed a molecular clock,
  // then half this tree distance would tell us how long ago the duplication
  // happened, and each of the two objects would be exactly this distance from
  // the duplication node.  Because the tree does not obey a molecular clock,
  // for any particular duplication node, the distances we get this way from
  // different pairs of duplicated objects that descend from it may be
  // different.  For each duplication node, we will record the maximum of all
  // such distances for any pair of duplicated objects descending from that
  // node, and which pair of objects yielded this maximum.  Call this maximum
  // duplication_node_distance.  

  // For each maximal node, it is maximal due to several possible pairs of
  // alternative nearest objects.  It may not always be the case that for each
  // pair of alternative nearest objects causing a particular node to be
  // maximal, one or both members of the pair will be descendants of the maximal
  // node.  The duplication node between them also may not be a descendant of
  // the maximal node.  The duplication distance of a maximal node will be the
  // maximum duplication_node_distance for all duplication nodes of pairs of
  // alternative nearest objects to the maximal node.  We will find the
  // duplication distance for each maximal node, and record which pair of
  // objects has the duplication_node_distance yielding that duplication
  // distance.  Note that this pair of objects may not be the same as the one
  // yielding the duplication_node_distance of the duplication node.

  // We find the distances from each of the objects that cause nodes to be
  // maximal to each of their ancestors.
  map<int, map<ObjectT, double> *> distances_to_maximizing_objects;
  find_distances_to_maximizing_objects(tree, breadth_first_visit_order,
                                    alternative_nearest_objects,
                                    leaf_of_object,
                                    distances_to_maximizing_objects);

  // The left_ids and right_ids of the nodes in this tree have previously been
  // computed and passed in idMap (which has idMap[nodeId] == pair<left_id,
  // right_id>).  These numbers, which were computed by a modified preorder tree
  // traversal, have the useful properties that:
  // node N1 descends from ndoe N2 if and only if (iff):
  //   left_id(N1) >= left_id(N2) and right_id(N1) <= right_id(N2)
  // node N is an ancestor of nodes N1, ..., Nk iff:
  //   left_id(Ni) >= left_id(N) for all i = 1,...,k and
  //   right_id(Ni) <= right_id(N) for all i = 1,...,k
  // node N is the most recent common ancestor of nodes N1, ..., Nk iff
  //   node N is an ancestor of nodes N1, ..., Nk and
  //   any of the following equivalent conditions holds:
  //   (i) node N has the greatest left_id among ancestors of nodes N1,...,Nk
  //   (ii) node N has the least right_id among ancestors of nodes N1,...,Nk
  //   (iii) node N appears last in the breadth_first_search_order among
  //   ancestors of nodes N1,...,Nk
  // We'll be using (iii).

  // We'll be looking up nodes by either their left_id or their right_id.  The
  // nodeIdOfLeftIdMap was passed in, but we now create a nodeIdOfRighgtIdMap as
  // well, using the idMap (which has idMap[nodeId] == pair<left_id, right_id>).
  map<int, int> nodeIdOfRightIdMap;
  map<int, pair<int, int> >::const_iterator node_id_iter;
  for (node_id_iter = idMap.begin(); node_id_iter != idMap.end();
        ++node_id_iter) {
    nodeIdOfRightIdMap[node_id_iter->second.second] = node_id_iter->first;
  }

  // We will be finding the tree distances between pairs of
  // alternative_nearest_objects.  To make it easy to test the condition that
  // N is an ancestor of N1 and N2, i.e., that
  // left_id(N1) >= left_id(N) and left_id(N2) >= left_id(N) and
  // right_id(N1) <= right_id(N) and right_id(N2) <= right_id(N),
  // we will keep track of each pair of leaves N1 and N2 as
  // left_id = min(left_id(N1), left_id(N2)) and
  // right_id = max(right_id(N1), right_id(N2)).
  map<int, vector<int> *> other_leaf_right_ids_of_leaf_left_ids;
  map<int, vector<int> *>::const_iterator right_left_iter;
  int leftId0, leftId1, rightId0, rightId1, leftId, rightId;
  map<ObjectT, set<int> *>::const_iterator obj_nodes_iter;
  map<ObjectT, map<ObjectT, set<int> *> *>::const_iterator
    obj_objs_nodes_iter;
  int obj0, obj1, lesserObj, greaterObj;
  map<ObjectT, int>::const_iterator obj_leaf_iter;
  map<int, pair<int,int> >::const_iterator node_left_right_iter;
  map<ObjectT, map<ObjectT, DupInfo *> *>::iterator obj_obj_dup_iter;
  for (obj_objs_nodes_iter 
        = alternative_nearest_objects.begin();
      obj_objs_nodes_iter
        != alternative_nearest_objects.end();
      ++obj_objs_nodes_iter) {
    obj0 = obj_objs_nodes_iter->first;
    obj_leaf_iter = leaf_of_object.find(obj0);
    node_left_right_iter = idMap.find(obj_leaf_iter->second);
    leftId0 = node_left_right_iter->second.first;
    rightId0 = node_left_right_iter->second.second;
    for (obj_nodes_iter = obj_objs_nodes_iter->second->begin();
        obj_nodes_iter != obj_objs_nodes_iter->second->end();
        ++obj_nodes_iter) {
      obj1 = obj_nodes_iter->first;
      lesserObj = (obj0 < obj1) ? obj0 : obj1;
      greaterObj = (obj0 < obj1) ? obj1 : obj0;
      obj_obj_dup_iter = duplication_node_of_objects.find(lesserObj);
      if (obj_obj_dup_iter == duplication_node_of_objects.end()) {
        duplication_node_of_objects[lesserObj] = new map<ObjectT, DupInfo *>;
      }
      obj_leaf_iter = leaf_of_object.find(obj1);
      node_left_right_iter = idMap.find(obj_leaf_iter->second);
      leftId1 = node_left_right_iter->second.first;
      rightId1 = node_left_right_iter->second.second;
      leftId = (leftId0 < leftId1) ? leftId0 : leftId1;
      rightId = (rightId0 > rightId1) ? rightId0 : rightId1;
      right_left_iter = other_leaf_right_ids_of_leaf_left_ids.find(leftId);
      if (right_left_iter == other_leaf_right_ids_of_leaf_left_ids.end()) {
        other_leaf_right_ids_of_leaf_left_ids[leftId] = new vector<int>;
      }
      other_leaf_right_ids_of_leaf_left_ids[leftId]->push_back(rightId);
    }
  }
  // We sort the left_ids and right_ids.  This way when we are checking whether
  // a node is an ancestor of some pair, it will be easy for us to only check
  // those pairs of which it could possibly be an ancestor.
  vector<int> sorted_left_ids_of_leaves;
  for (right_left_iter = other_leaf_right_ids_of_leaf_left_ids.begin(); 
      right_left_iter != other_leaf_right_ids_of_leaf_left_ids.end(); 
      ++right_left_iter) {
    sorted_left_ids_of_leaves.push_back(right_left_iter->first);
    std::sort(right_left_iter->second->begin(),
              right_left_iter->second->end());
    std::reverse(right_left_iter->second->begin(),
              right_left_iter->second->end());
  }
  std::sort(sorted_left_ids_of_leaves.begin(), sorted_left_ids_of_leaves.end());

  int least_left_id_greater_than_node_left_id,
      greatest_right_id_less_than_node_right_id;
  vector<int>::iterator j, k;
  int r;
  AttributeT attr;
  int nodeId, childId;
  double distanceToLesserObj, distanceToGreaterObj, duplication_node_distance,
          new_duplication_node_distance;
  map<int, ObjObjInfo *>::iterator node_obj_pair_iter;
  map<int, int>::const_iterator int_node_iter;
  map<int, AttributeT>::const_iterator node_attr_iter;
  map<int, ObjectT>::const_iterator leaf_obj_iter;
  map<ObjectT, int>::const_iterator ancestral_iter;
  map<int, map<AttributeT, TreeDistanceInfo<ObjectT,
                  nullObjectValue> *> *>::const_iterator dist_attr_node_iter;
  // Now we go up the tree computing the distances between each pair of
  // alternative nearest objects, and all the duplication_node_distances.
  for (int i = breadth_first_visit_order.size() - 1; i >= 0; --i) {
    nodeId = breadth_first_visit_order[i];
    const vector<int> &children = tree.getSonsId(nodeId);
    duplication_node_distance = -1.0;
    for (int c = 0; c < children.size(); ++c) {
      childId = children[c];
      node_obj_pair_iter
        = pair_of_objects_yielding_duplication_node_distance.find(childId);
      if (node_obj_pair_iter
          != pair_of_objects_yielding_duplication_node_distance.end()) {
        new_duplication_node_distance = getDuplicationNodeDistance(childId,
                          duplication_node_of_objects,
                        pair_of_objects_yielding_duplication_node_distance);
        if (duplication_node_distance < 0.0) {
          pair_of_objects_yielding_duplication_node_distance[nodeId]
            = new ObjObjInfo();
          pair_of_objects_yielding_duplication_node_distance[
                              nodeId]->setValues(
            pair_of_objects_yielding_duplication_node_distance[
                                                    childId]->getLesserObj(),
            pair_of_objects_yielding_duplication_node_distance[
                                                    childId]->getGreaterObj());
          duplication_node_distance = new_duplication_node_distance;
        } else {
          if (new_duplication_node_distance > duplication_node_distance) {
            pair_of_objects_yielding_duplication_node_distance[
                                nodeId]->setValues(
              pair_of_objects_yielding_duplication_node_distance[
                                                    childId]->getLesserObj(),
              pair_of_objects_yielding_duplication_node_distance[
                                                    childId]->getGreaterObj());
            duplication_node_distance = new_duplication_node_distance;
          }
        }
      }
    }
    node_left_right_iter = idMap.find(nodeId);
    leftId = node_left_right_iter->second.first;
    rightId = node_left_right_iter->second.second;
    // If the binary_search returns -1, this means leftId is before the
    // beginning of sorted_left_ids_of_leaves.  So it will be correct that j =
    // sorted_left_ids_of_leaves.begin()
    j = sorted_left_ids_of_leaves.begin() + 1 
            + binary_search(leftId, sorted_left_ids_of_leaves);
    least_left_id_greater_than_node_left_id = *j;
    // least_left_id_greater_than_node_left_id will be increasing as we go
    // through this loop; we stop when it goes past rightId.
    while (j != sorted_left_ids_of_leaves.end()
          && least_left_id_greater_than_node_left_id < rightId) {
      int_node_iter 
        = nodeIdOfLeftIdMap.find(least_left_id_greater_than_node_left_id);
      node_attr_iter = attribute_of_node.find(int_node_iter->second);
      leaf_obj_iter = object_of_leaf.find(int_node_iter->second);
      obj0 = leaf_obj_iter->second;
      // If the binary_search returns -1, this means rightId is after the end of
      // the right_ids that are paired with
      // least_left_id_greater_than_node_left_id (since they are sorted in
      // reverse order).  So it will be correct that k =
      // other_leaf_right_ids_of_leaf_left_ids[least_left_id_greater_than_node_left_id]->begin().
      // beginning of the right_ids that are paired with
      // least_left_id_greater_than_node_left_id.  Hence none of them can be
      // descendant from this node, and we should continue.
      k = other_leaf_right_ids_of_leaf_left_ids[
                              least_left_id_greater_than_node_left_id]->begin()
          + 1 + reverse_binary_search(rightId,
                  *(other_leaf_right_ids_of_leaf_left_ids[
                                  least_left_id_greater_than_node_left_id]));
      greatest_right_id_less_than_node_right_id = *k;
      // greatest_right_id_less_than_node_right_id will be decreasing as we go
      // through this loop; we stop when it goes past leftId.
      while (k != other_leaf_right_ids_of_leaf_left_ids[
                          least_left_id_greater_than_node_left_id]->end()
          && greatest_right_id_less_than_node_right_id > leftId) {
        int_node_iter = nodeIdOfRightIdMap.find(
                                  greatest_right_id_less_than_node_right_id);
        leaf_obj_iter = object_of_leaf.find(int_node_iter->second);
        obj1 = leaf_obj_iter->second;
        // This node is the duplication node of obj0, obj1
        lesserObj = (obj0 < obj1) ? obj0 : obj1;
        greaterObj = (obj0 < obj1) ? obj1 : obj0;
        distanceToLesserObj 
          = (*distances_to_maximizing_objects[nodeId])[lesserObj];
        distanceToGreaterObj
          = (*distances_to_maximizing_objects[nodeId])[greaterObj];
        (*duplication_node_of_objects[lesserObj])[greaterObj]
          = new DupInfo(nodeId, distanceToLesserObj, distanceToGreaterObj);
        new_duplication_node_distance 
          = (distanceToLesserObj + distanceToGreaterObj) / 2;
        if (duplication_node_distance < 0.0) {
          pair_of_objects_yielding_duplication_node_distance[nodeId]
            = new ObjObjInfo();
          pair_of_objects_yielding_duplication_node_distance[
                              nodeId]->setValues(lesserObj, greaterObj);
          duplication_node_distance = new_duplication_node_distance;
        } else {
          if (new_duplication_node_distance > duplication_node_distance) {
            pair_of_objects_yielding_duplication_node_distance[
                                nodeId]->setValues(lesserObj, greaterObj);
            duplication_node_distance = new_duplication_node_distance;
          }
        }
        // Since we've found the MRCA of obj0, obj1, we don't need to check
        // for this pair any more.  In fact, we had better not, otherwise we
        // would erroneously keep setting the duplication node to other
        // ancestors closer to the root of the tree.  So we must erase
        // right_id here.
        k = other_leaf_right_ids_of_leaf_left_ids[
                          least_left_id_greater_than_node_left_id]->erase(k);
        if (k != other_leaf_right_ids_of_leaf_left_ids[
                            least_left_id_greater_than_node_left_id]->end()) {
          // Now k points to the item that originally came *after*
          // greatest_right_id_less_than_node_right_id in the vector, which
          // means that it is *less* than
          // greatest_right_id_less_than_node_right_id (since the vector is in
          // reverse order), and hence still less than rightId.
          greatest_right_id_less_than_node_right_id = *k;
        }
      }
      if (other_leaf_right_ids_of_leaf_left_ids[
                    least_left_id_greater_than_node_left_id]->size() == 0) {
        // Since we've found the duplication nodes of all pairs with obj0, we
        // don't need to check it any more.  So we can erase the left_id.
        delete other_leaf_right_ids_of_leaf_left_ids[
                    least_left_id_greater_than_node_left_id];
        other_leaf_right_ids_of_leaf_left_ids.erase(
                                      least_left_id_greater_than_node_left_id);
        j = sorted_left_ids_of_leaves.erase(j);
      } else {
        ++j;
      }
      // Now j points to the item that originall came after
      // least_left_id_greater_than_node_left_id in sorted_left_ids_of_leaves,
      // which means that it greater than
      // least_left_id_greater_than_node_left_id, and hence also than leftId.
      if (j != sorted_left_ids_of_leaves.end()) {
        least_left_id_greater_than_node_left_id = *j;
      }
    }
  }
  if (other_leaf_right_ids_of_leaf_left_ids.size() > 0) {
    cerr << "Warning: other_leaf_right_ids_of_leaf_left_ids ended nonempty"
      << endl;
  }

  // Now we go through alternative_nearest_objects computing the duplication
  // distances of the maximal nodes.
  set<int>::iterator maximal_node_iter;
  int duplicationNodeId;
  map<int, ObjObjInfo *>::iterator obj_obj_max_node_iter;
  for (obj_objs_nodes_iter 
        = alternative_nearest_objects.begin();
      obj_objs_nodes_iter
        != alternative_nearest_objects.end();
      ++obj_objs_nodes_iter) {
    obj0 = obj_objs_nodes_iter->first;
    for (obj_nodes_iter = obj_objs_nodes_iter->second->begin();
        obj_nodes_iter != obj_objs_nodes_iter->second->end();
        ++obj_nodes_iter) {
      obj1 = obj_nodes_iter->first;
      lesserObj = (obj0 < obj1) ? obj0 : obj1;
      greaterObj = (obj0 < obj1) ? obj1 : obj0;
      duplicationNodeId 
        = (*duplication_node_of_objects[lesserObj])
                                          [greaterObj]->getDuplicationNodeId();
      for (maximal_node_iter = obj_nodes_iter->second->begin();
            maximal_node_iter != obj_nodes_iter->second->end();
            ++maximal_node_iter) {
        obj_obj_max_node_iter = pair_of_nearest_objects_to_maximal_node.find(
                                                          *maximal_node_iter);
        if (obj_obj_max_node_iter 
              == pair_of_nearest_objects_to_maximal_node.end()) {
          pair_of_nearest_objects_to_maximal_node[*maximal_node_iter]
            = new ObjObjInfo();
          pair_of_nearest_objects_to_maximal_node[
                        *maximal_node_iter]->setValues(lesserObj, greaterObj);

        } else {
          if (getDuplicationNodeDistance(duplicationNodeId,
                            duplication_node_of_objects,
                            pair_of_objects_yielding_duplication_node_distance)
              > getMaximalNodeDistance(*maximal_node_iter,
                            duplication_node_of_objects,
                            pair_of_objects_yielding_duplication_node_distance,
                            pair_of_nearest_objects_to_maximal_node)) {
            pair_of_nearest_objects_to_maximal_node[
                        *maximal_node_iter]->setValues(lesserObj, greaterObj);
          }
        }
      }
    }
  }

  // Now we go up the tree finding, for each maximal node, the greatest
  // duplication distance of any of its strict descendants (i.e., not including
  // itself) which is also a maximal node.  We store the pair of the maximal
  // node's own duplication distance and this for each of the maximal nodes, and
  // for the root.
  map<int, double> all_node_greatest_distance_of_maximal_descendant;
  map<int, double>::const_iterator node_iter, parent_iter;
  int parentId;
  double currentDistance;
  for (int i = breadth_first_visit_order.size() - 1; i >= 1; --i) {
    nodeId = breadth_first_visit_order[i];
    node_iter = all_node_greatest_distance_of_maximal_descendant.find(nodeId);
    parentId = tree.getFatherId(nodeId);
    parent_iter 
      = all_node_greatest_distance_of_maximal_descendant.find(parentId);
    if (node_iter != all_node_greatest_distance_of_maximal_descendant.end()) {
      if (parent_iter 
          == all_node_greatest_distance_of_maximal_descendant.end()) {
        all_node_greatest_distance_of_maximal_descendant[parentId]
          = node_iter->second;
      } else {
        if (node_iter->second > parent_iter->second) {
          all_node_greatest_distance_of_maximal_descendant[parentId]
            = node_iter->second;
        }
      }
    }
    obj_obj_max_node_iter 
      = pair_of_nearest_objects_to_maximal_node.find(nodeId);
    if (obj_obj_max_node_iter 
        != pair_of_nearest_objects_to_maximal_node.end()) {
      if (node_iter == all_node_greatest_distance_of_maximal_descendant.end()) {
        greatest_distance_of_maximal_descendant[nodeId] = -1.0;
      } else {
        greatest_distance_of_maximal_descendant[nodeId] = node_iter->second;
      }
      currentDistance = getMaximalNodeDistance(nodeId, 
                            duplication_node_of_objects,
                            pair_of_objects_yielding_duplication_node_distance,
                            pair_of_nearest_objects_to_maximal_node);
      if (parent_iter 
          == all_node_greatest_distance_of_maximal_descendant.end()) {
        all_node_greatest_distance_of_maximal_descendant[parentId]
          = currentDistance;
      } else {
        if (currentDistance > parent_iter->second) {
          all_node_greatest_distance_of_maximal_descendant[parentId]
            = currentDistance;
        }
      }
    }
  }
  int rootId = tree.getRootId();
  obj_obj_max_node_iter
    = pair_of_nearest_objects_to_maximal_node.find(rootId);
  node_iter = all_node_greatest_distance_of_maximal_descendant.find(rootId);
  if (node_iter == all_node_greatest_distance_of_maximal_descendant.end()) {
    greatest_distance_of_maximal_descendant[rootId] = -1.0;
  } else {
    greatest_distance_of_maximal_descendant[rootId] = node_iter->second;
  }
      
  map<int, map<ObjectT, double> *>::iterator node_obj_dist_iter; 
  for (node_obj_dist_iter = distances_to_maximizing_objects.begin();
      node_obj_dist_iter != distances_to_maximizing_objects.end();
      ++node_obj_dist_iter) {
    delete node_obj_dist_iter->second;
  }
}
