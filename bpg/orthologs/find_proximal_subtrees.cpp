// find_proximal_subtrees.cpp
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

#include "modified_preorder_tree_traversal.h"
#include "find_proximal_subtrees.h"
#include <iostream>

//template<class ObjectT, ObjectT nullObjectValue, class AttributeT>
void addMaximalNode(
        const LeftRightIdsOfNodeIdMap &idMap,
        const map<ObjectT, int> leaf_of_object,
        const set<ObjectT> &nearest_objects_to_maximal_node_neighbors,
        map<int, set<AttributeT > *> &attributes_for_which_node_is_maximal,
        map<ObjectT , set<int> *> &maximal_nodes_of_object,
        map<ObjectT, map<ObjectT, set<int> *> *> &alternative_nearest_objects,
        int nodeId, AttributeT attr, ObjectT obj) {
//  typename 
  map<int, set<AttributeT > *>::const_iterator attrs_node_iter;
//  typename 
  map<AttributeT , set<int> *>::const_iterator nodes_attr_iter;
//  typename 
  map<ObjectT , set<int> *>::const_iterator nodes_obj_iter;
  attrs_node_iter = attributes_for_which_node_is_maximal.find(nodeId);
  if (attrs_node_iter == attributes_for_which_node_is_maximal.end()) {
    attributes_for_which_node_is_maximal[nodeId] = new set<AttributeT >;
  }
  attributes_for_which_node_is_maximal[nodeId]->insert(attr);
  nodes_obj_iter = maximal_nodes_of_object.find(obj);
  if (nodes_obj_iter == maximal_nodes_of_object.end()) {
    maximal_nodes_of_object[obj] = new set<int>;
  }
  maximal_nodes_of_object[obj]->insert(nodeId);
  map<ObjectT, int>::const_iterator leaf_obj_iter;
  leaf_obj_iter = leaf_of_object.find(obj);
  int leafId = (*leaf_obj_iter).second;
  LeftRightIdsOfNodeIdMap::const_iterator leaf_left_iter, node_left_iter;
  int leaf_left_id, leaf_right_id, node_left_id, node_right_id;
  node_left_iter = idMap.find(nodeId);
  node_left_id = ((*node_left_iter).second).first;
  node_right_id = ((*node_left_iter).second).second;
  leaf_left_iter = idMap.find(leafId);
  leaf_left_id = ((*leaf_left_iter).second).first;
  leaf_right_id = ((*leaf_left_iter).second).second;
  set<ObjectT>::const_iterator alt_objs_iter;
  map<ObjectT, map<ObjectT, set<int> *> *>::iterator obj_objs_nodes_iter;
  map<ObjectT, set<int> *>::iterator obj_nodes_iter;
  ObjectT obj0, obj1;
  for (alt_objs_iter = nearest_objects_to_maximal_node_neighbors.begin();
      alt_objs_iter != nearest_objects_to_maximal_node_neighbors.end();
      ++alt_objs_iter) {
    if (obj == *alt_objs_iter) {
      continue;
    }
    if (obj < *alt_objs_iter) {
      obj0 = obj;
      obj1 = *alt_objs_iter;
    } else {
      obj0 = *alt_objs_iter;
      obj1 = obj;
    }
    obj_objs_nodes_iter = alternative_nearest_objects.find(obj0);
    if (obj_objs_nodes_iter == alternative_nearest_objects.end()) {
      alternative_nearest_objects[obj0] = new map<ObjectT, set<int> *>;
    }
    obj_nodes_iter = alternative_nearest_objects[obj0]->find(obj1);
    if (obj_nodes_iter == alternative_nearest_objects[obj0]->end()) {
      (*alternative_nearest_objects[obj0])[obj1] = new set<int>;
    }
    (*alternative_nearest_objects[obj0])[obj1]->insert(nodeId);
  }
}

//template<class ObjectT, ObjectT nullObjectValue, class AttributeT>
void find_proximal_subtrees(const bpp::Tree &tree,
        const vector<int> &breadth_first_visit_order,
        const map<int, AttributeT> attribute_of_node,
        const map<AttributeT , ObjectT> unique_object_with_attribute,
        const set<AttributeT> &attributes_with_multiple_objects,
        const map<int, map<AttributeT, TreeDistanceInfo<ObjectT, 
                                                    nullObjectValue> *> *>
          distance_from_object_with_attribute_to_node,
        const LeftRightIdsOfNodeIdMap &idMap,
        const map<ObjectT, int> &leaf_of_object,
        const map<int, ObjectT> &object_of_leaf,
        map<int, set<AttributeT> *> 
          &attributes_with_multiple_objects_for_which_node_is_maximal,
        set<int> &super_orthologous_nodes,
        map<ObjectT , set<int> *> &maximal_nodes_of_object,
        map<ObjectT, map<ObjectT, set<int> *> *> &alternative_nearest_objects) {
  // DO NOT clear attribute_of_node, unique_object_with_attribute, or
  // distance_from_object_with_attribute_to_node.  These have already been
  // filled, except for the isObjectPure field of the TreeDistanceInfo's.
  attributes_with_multiple_objects_for_which_node_is_maximal.clear();
  super_orthologous_nodes.clear();
  maximal_nodes_of_object.clear();
  map<int, bool> has_noninparalogous_maximal_descendant;
//  typename 
  map<int, map<AttributeT, TreeDistanceInfo<ObjectT,
                  nullObjectValue> *> *>::const_iterator dist_attr_node_iter,
                                                    child_dist_attr_node_iter;
//  typename 
  map<AttributeT, 
      TreeDistanceInfo<ObjectT, nullObjectValue> *>::const_iterator 
        attr_node_iter, child_attr_node_iter;
  int rootId = tree.getRootId();
//  typename 
  set<AttributeT >::const_iterator attr_iter;
//  typename 
  map<AttributeT , ObjectT >::const_iterator obj_attr_iter;
//  typename 
  map<ObjectT , set<int> *>::const_iterator nodes_obj_iter;
  for (obj_attr_iter = unique_object_with_attribute.begin();
        obj_attr_iter != unique_object_with_attribute.end();
        ++obj_attr_iter) {
    const AttributeT attr = (*obj_attr_iter).first;
    const ObjectT obj = (*obj_attr_iter).second;
    nodes_obj_iter = maximal_nodes_of_object.find(obj);
    if (nodes_obj_iter == maximal_nodes_of_object.end()) {
      maximal_nodes_of_object[obj] = new set<int>;
    }
    maximal_nodes_of_object[obj]->insert(rootId);
  }
  if (attributes_with_multiple_objects.empty() 
      || tree.getNumberOfNodes() <= 1) {
    if (attribute_of_node.find(rootId) == attribute_of_node.end()
        && unique_object_with_attribute.size() > 0) {
      cout << "Root is not inparalogous and is only maximal node "
            << "=> superorthologous" << endl;
      super_orthologous_nodes.insert(rootId);
    }
    return;
  }
  int nodeId, childId;
  bool isObjectPure, isChildObjectPure;
  map<int, map<AttributeT, set<ObjectT> *> *> 
    nearest_objects_to_node_neighbors;
  map<int, AttributeT>::const_iterator attr_leaf_iter;
  for(int i = breadth_first_visit_order.size() - 1; i >= 0; --i) {
    nodeId = breadth_first_visit_order[i];
    dist_attr_node_iter 
      = distance_from_object_with_attribute_to_node.find(nodeId);
    has_noninparalogous_maximal_descendant[nodeId] = false;
    nearest_objects_to_node_neighbors[nodeId] 
      = new map<AttributeT, set<ObjectT> *>;
    if (tree.isLeaf(nodeId)) {
      for(attr_iter = attributes_with_multiple_objects.begin();
          attr_iter != attributes_with_multiple_objects.end();
          ++attr_iter) {
        attr_node_iter = ((*dist_attr_node_iter).second)->find(*attr_iter);
        ((*attr_node_iter).second)->setIsObjectPure(true);
        (*nearest_objects_to_node_neighbors[nodeId])[*attr_iter] 
          = new set<ObjectT>;
      }
      attr_leaf_iter = attribute_of_node.find(nodeId);
      if (attr_leaf_iter != attribute_of_node.end() 
          && attributes_with_multiple_objects.find(attr_leaf_iter->second)
              != attributes_with_multiple_objects.end()) {
        (*nearest_objects_to_node_neighbors[nodeId])
          [attr_leaf_iter->second]->insert(
            object_of_leaf.find(nodeId)->second);
      }
      continue;
    }
    const vector<int> &children = tree.getSonsId(nodeId);
    if (children.size() > 0) {
      for (int j = 0; j < children.size(); ++j) {
        has_noninparalogous_maximal_descendant[nodeId] 
            = has_noninparalogous_maximal_descendant[nodeId]
              || has_noninparalogous_maximal_descendant[children[j]];
      }
      for(attr_iter = attributes_with_multiple_objects.begin();
          attr_iter != attributes_with_multiple_objects.end();
          ++attr_iter) {
        isObjectPure = true;
        (*nearest_objects_to_node_neighbors[nodeId])[*attr_iter]
          = new set<ObjectT>;
        attr_node_iter = ((*dist_attr_node_iter).second)->find(*attr_iter);
        const ObjectT nearest_object
          = ((*attr_node_iter).second)->getNearestObjectWithAttributeValue();
        (*nearest_objects_to_node_neighbors[nodeId])
                                        [*attr_iter]->insert(nearest_object);
        for (int j = 0; j < children.size(); ++j) {
          childId = children[j];
          child_dist_attr_node_iter 
            = distance_from_object_with_attribute_to_node.find(childId);
          child_attr_node_iter 
            = ((*child_dist_attr_node_iter).second)->find(*attr_iter);
          const ObjectT nearest_object_to_child 
            = ((*child_attr_node_iter).second)
                ->getNearestObjectWithAttributeValue();
          if (attribute_of_node.find(childId) == attribute_of_node.end() ||
              (*(attribute_of_node.find(childId))).second != *attr_iter) {
            // The child node is not inparalogous, or if it is, it's for a
            // different species than we're considering now.
            for (set<ObjectT>::const_iterator child_obj_iter = 
                  (*nearest_objects_to_node_neighbors[childId])
                  [*attr_iter]->begin();
                  child_obj_iter != (*nearest_objects_to_node_neighbors[childId])
                    [*attr_iter]->end();
                  ++child_obj_iter) {
                (*nearest_objects_to_node_neighbors[nodeId])
                                          [*attr_iter]->insert(*child_obj_iter);
            }
          } else {
            // The child node is inparalogous with the same species we're
            // considering now.  So only consider the nearest sequence to the
            // child itself (not any of the duplicated nearest sequences
            // underneath it).
            (*nearest_objects_to_node_neighbors[nodeId])
                                [*attr_iter]->insert(nearest_object_to_child);
          }
          isChildObjectPure 
            = ((*child_attr_node_iter).second)->isObjectPure();
          if (!isChildObjectPure || 
              nearest_object_to_child != nearest_object &&
              attribute_of_node.find(nodeId) == attribute_of_node.end()) {
            ((*attr_node_iter).second)->setIsObjectPure(false);
            isObjectPure = false;
          }
        }
        ((*attr_node_iter).second)->setIsObjectPure(isObjectPure);
        if (!isObjectPure) {
          for (int j = 0; j < children.size(); ++j) {
            childId = children[j];
            child_dist_attr_node_iter 
              = distance_from_object_with_attribute_to_node.find(childId);
            child_attr_node_iter 
              = ((*child_dist_attr_node_iter).second)->find(*attr_iter);
            const ObjectT obj 
              = ((*child_attr_node_iter).second)
                  ->getNearestObjectWithAttributeValue();
            isChildObjectPure 
              = ((*child_attr_node_iter).second)->isObjectPure();
            if (isChildObjectPure) {
              addMaximalNode(
                  idMap, leaf_of_object,
                  *((*nearest_objects_to_node_neighbors[nodeId])[*attr_iter]),
                  attributes_with_multiple_objects_for_which_node_is_maximal,
                            maximal_nodes_of_object, 
                            alternative_nearest_objects,
                            childId, *attr_iter, obj);
              if (attribute_of_node.find(childId) 
                  == attribute_of_node.end()) {
                if (!has_noninparalogous_maximal_descendant[childId]) {
                  super_orthologous_nodes.insert(childId);
                }
                has_noninparalogous_maximal_descendant[nodeId] = true;
              } else if ((*(attribute_of_node.find(childId))).second 
                          != *attr_iter) {
                has_noninparalogous_maximal_descendant[nodeId] = true;
              }
            }
          }
        }
      }
    }
  }
  dist_attr_node_iter 
    = distance_from_object_with_attribute_to_node.find(rootId);
  for(attr_iter = attributes_with_multiple_objects.begin();
      attr_iter != attributes_with_multiple_objects.end();
      ++attr_iter) {
    attr_node_iter = ((*dist_attr_node_iter).second)->find(*attr_iter);
    const ObjectT obj
      = ((*attr_node_iter).second)->getNearestObjectWithAttributeValue();
    isObjectPure = ((*attr_node_iter).second)->isObjectPure();
    if (isObjectPure) {
      addMaximalNode(
          idMap, leaf_of_object,
          // This should only have one object
          // The only way that the root can be object pure for an attribute with
          // multiple objects in the tree is if all the objects for that
          // attribute occur in an inparalogous subtree.  So this just records
          // which of those inparalogous objects has the shortest distance from
          // the root.
          *((*nearest_objects_to_node_neighbors[rootId])[*attr_iter]),
          attributes_with_multiple_objects_for_which_node_is_maximal,
                    maximal_nodes_of_object, 
                    alternative_nearest_objects, rootId, *attr_iter, obj);
    }
  }
  if (!has_noninparalogous_maximal_descendant[rootId] 
      && attribute_of_node.find(rootId) == attribute_of_node.end()
      && (attributes_with_multiple_objects_for_which_node_is_maximal.find(
                                                                      rootId) 
        != attributes_with_multiple_objects_for_which_node_is_maximal.end()
        || unique_object_with_attribute.size() > 0)) {
    super_orthologous_nodes.insert(rootId);
  }
}
