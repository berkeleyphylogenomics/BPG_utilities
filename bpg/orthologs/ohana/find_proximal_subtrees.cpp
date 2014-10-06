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

#include "find_proximal_subtrees.h"
#include <iostream>

//template<class ObjectT, ObjectT nullObjectValue, class AttributeT>
void addMaximalNode(
        map<int, set<AttributeT > *> &attributes_for_which_node_is_maximal,
        map<ObjectT , set<int> *> &maximal_nodes_of_object,
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
        map<int, set<AttributeT> *> 
          &attributes_with_multiple_objects_for_which_node_is_maximal,
        set<int> &super_orthologous_nodes,
        map<ObjectT , set<int> *> &maximal_nodes_of_object) {
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
  for(int i = breadth_first_visit_order.size() - 1; i >= 0; --i) {
    nodeId = breadth_first_visit_order[i];
    dist_attr_node_iter 
      = distance_from_object_with_attribute_to_node.find(nodeId);
    has_noninparalogous_maximal_descendant[nodeId] = false;
    if (tree.isLeaf(nodeId)) {
      for(attr_iter = attributes_with_multiple_objects.begin();
          attr_iter != attributes_with_multiple_objects.end();
          ++attr_iter) {
        attr_node_iter = ((*dist_attr_node_iter).second)->find(*attr_iter);
        ((*attr_node_iter).second)->setIsObjectPure(true);
      }
      continue;
    }
    vector<int> children = tree.getSonsId(nodeId);
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
        attr_node_iter = ((*dist_attr_node_iter).second)->find(*attr_iter);
        const ObjectT nearest_object
          = ((*attr_node_iter).second)->getNearestObjectWithAttributeValue();
        for (int j = 0; j < children.size() && isObjectPure; ++j) {
          childId = children[j];
          child_dist_attr_node_iter 
            = distance_from_object_with_attribute_to_node.find(childId);
          child_attr_node_iter 
            = ((*child_dist_attr_node_iter).second)->find(*attr_iter);
          const ObjectT nearest_object_to_child 
            = ((*child_attr_node_iter).second)
                ->getNearestObjectWithAttributeValue();
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
                  attributes_with_multiple_objects_for_which_node_is_maximal,
                            maximal_nodes_of_object, 
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
          attributes_with_multiple_objects_for_which_node_is_maximal,
                    maximal_nodes_of_object, rootId, *attr_iter, obj);
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
