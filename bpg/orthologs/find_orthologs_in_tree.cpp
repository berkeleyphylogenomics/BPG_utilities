// find_orthologs_in_tree.h
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

#include "find_orthologs_in_tree.h"

using namespace std;
using namespace bpp;

//template<class ObjectT, ObjectT nullObjectValue, class AttributeT, 
//          AttributeT nullAttributeValue>
void find_orthologs_in_tree(const bpp::Tree &tree, 
        const vector<int> &breadth_first_visit_order,
        const map<int, ObjectT > &object_of_leaf,
        const map<int, AttributeT> &attribute_of_node,
        map<AttributeT, ObjectT > &unique_object_with_attribute,
        set<AttributeT> &attributes_with_multiple_objects,
        map<int, map<AttributeT, TreeDistanceInfo<ObjectT, 
                                                    nullObjectValue> *> *>
          &distance_from_object_with_attribute_to_node,
        bool suppress_inparalogs, bool seek_longest) {
  distance_from_object_with_attribute_to_node.clear();
  unique_object_with_attribute.clear();
  const vector<int> &leafIds = tree.getLeavesId();
  const vector<int> &nodeIds = tree.getNodesId();
//  typename 
  map<AttributeT, ObjectT>::const_iterator 
      obj_attr_iter;
  vector<int>::const_iterator node_iter;
  typedef map<int, AttributeT > AttributeNodeMap;
//  typename 
  AttributeNodeMap::const_iterator attribute_of_node_iter;
  typedef map<int, ObjectT > ObjectLeafMap;
//  typename
  ObjectLeafMap::const_iterator object_of_leaf_iter;
  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    attribute_of_node_iter = attribute_of_node.find(*node_iter);
    const AttributeT attributeValue = (*attribute_of_node_iter).second;
    object_of_leaf_iter = object_of_leaf.find(*node_iter);
    const ObjectT object = (*object_of_leaf_iter).second;
    if (attributeValue != nullAttributeValue) {
      obj_attr_iter = unique_object_with_attribute.find(attributeValue);
      if (obj_attr_iter == unique_object_with_attribute.end() 
          && object != nullObjectValue) {
        unique_object_with_attribute[attributeValue] = object;
      } else {
        unique_object_with_attribute[attributeValue] = nullObjectValue;
        attributes_with_multiple_objects.insert(attributeValue);
      }
    }
  }
//  typename 
  set<AttributeT >::const_iterator attr_iter;
  for (attr_iter = attributes_with_multiple_objects.begin();
        attr_iter != attributes_with_multiple_objects.end();
        ++attr_iter) {
    unique_object_with_attribute.erase(*attr_iter);
  }
  if (attributes_with_multiple_objects.empty() 
      || tree.getNumberOfNodes() <= 1) {
    return;
  }

  for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
    distance_from_object_with_attribute_to_node[*node_iter]
      = new map<AttributeT , TreeDistanceInfo<ObjectT, 
                                                nullObjectValue> *>;
    for (attr_iter = attributes_with_multiple_objects.begin();
        attr_iter != attributes_with_multiple_objects.end();
        ++attr_iter) {
      (*(distance_from_object_with_attribute_to_node[*node_iter]))
        [*attr_iter] = new TreeDistanceInfo<ObjectT, nullObjectValue>();
    }
  }

  for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
    attribute_of_node_iter = attribute_of_node.find(*node_iter);
    const AttributeT attr = (*attribute_of_node_iter).second;
    obj_attr_iter = unique_object_with_attribute.find(attr);
    object_of_leaf_iter = object_of_leaf.find(*node_iter);
    if (attr != nullAttributeValue &&
        obj_attr_iter == unique_object_with_attribute.end() &&
        (*object_of_leaf_iter).second != nullObjectValue) {
      TreeDistanceInfo<ObjectT, nullObjectValue> &treeDistanceInfo = 
        *((*(distance_from_object_with_attribute_to_node[*node_iter]))[attr]);
      treeDistanceInfo.setNearestObjectWithAttributeValue(
                                              (*object_of_leaf_iter).second);
      treeDistanceInfo.setDistanceToNearestObjectWithAttributeValue(0.0);
      treeDistanceInfo.setNearestObjectIsDescendant(true);
    }
  }

  double attribute_to_node_distance, 
          attribute_to_parent_distance_through_node,
          current_attribute_to_parent_distance,
          attribute_to_node_next_nearest_distance,
          attribute_to_parent_next_nearest_distance_through_node,
          current_attribute_to_parent_next_nearest_distance;
  int nodeId, parentId;
  bool have_updated_next_nearest_from_node_nearest;
  for (int i = breadth_first_visit_order.size() - 1; i >= 1; --i) {
    nodeId = breadth_first_visit_order[i];
    parentId = tree.getFatherId(nodeId);
    attribute_of_node_iter = attribute_of_node.find(nodeId);
    for (attr_iter = attributes_with_multiple_objects.begin();
        attr_iter != attributes_with_multiple_objects.end(); 
        ++attr_iter) {
      have_updated_next_nearest_from_node_nearest = false;
      attribute_to_node_distance
        = (*(distance_from_object_with_attribute_to_node[nodeId]))
          [*attr_iter]->getDistanceToNearestObjectWithAttributeValue();
      if (attribute_to_node_distance != -1.0) {
        attribute_to_parent_distance_through_node 
          = attribute_to_node_distance + tree.getDistanceToFather(nodeId);
        current_attribute_to_parent_distance
          = (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->getDistanceToNearestObjectWithAttributeValue();
        const ObjectT nearest_object
          = (*(distance_from_object_with_attribute_to_node[nodeId]))
            [*attr_iter]->getNearestObjectWithAttributeValue();
        ObjectT current_nearest_object_to_parent
          = (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->getNearestObjectWithAttributeValue();
        const ObjectT current_next_nearest_object_to_parent
          = (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->getNextNearestObjectWithAttributeValue();
        current_attribute_to_parent_next_nearest_distance
          = (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->getDistanceToNextNearestObjectWithAttributeValue();
        if (current_attribute_to_parent_distance == -1.0 ||
              (seek_longest 
              ? (attribute_to_parent_distance_through_node
                      > current_attribute_to_parent_distance)
              : (attribute_to_parent_distance_through_node 
                      < current_attribute_to_parent_distance)) ||
              attribute_to_parent_distance_through_node
                == current_attribute_to_parent_distance
              && (current_nearest_object_to_parent == nullObjectValue
                  || nearest_object < current_nearest_object_to_parent)) {
          if (current_attribute_to_parent_distance != -1.0 &&
              current_nearest_object_to_parent != nullObjectValue &&
              current_nearest_object_to_parent != nearest_object &&
              (!suppress_inparalogs ||
              attribute_of_node_iter == attribute_of_node.end() 
              || (*attribute_of_node_iter).second != *attr_iter)) {
            (*(distance_from_object_with_attribute_to_node[parentId]))
              [*attr_iter]->setNextNearestObjectWithAttributeValue(
                            current_nearest_object_to_parent);
            (*(distance_from_object_with_attribute_to_node[parentId]))
              [*attr_iter]->setDistanceToNextNearestObjectWithAttributeValue(
                            current_attribute_to_parent_distance);
            current_attribute_to_parent_next_nearest_distance
              = current_attribute_to_parent_distance;
          }
          (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->setNearestObjectWithAttributeValue(nearest_object);
          current_nearest_object_to_parent = nearest_object;
          (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->setDistanceToNearestObjectWithAttributeValue(
              attribute_to_parent_distance_through_node);
          (*(distance_from_object_with_attribute_to_node[parentId]))
            [*attr_iter]->setNearestObjectIsDescendant(true);
        } else {
          if (nearest_object != current_nearest_object_to_parent &&
              (current_attribute_to_parent_next_nearest_distance == -1.0 ||
              ((seek_longest
                ? (attribute_to_parent_distance_through_node
                    > current_attribute_to_parent_next_nearest_distance)
                : (attribute_to_parent_distance_through_node
                    < current_attribute_to_parent_next_nearest_distance)) ||
               attribute_to_parent_distance_through_node
                == current_attribute_to_parent_next_nearest_distance &&
               nearest_object < current_next_nearest_object_to_parent) &&
              (!suppress_inparalogs ||
              attribute_of_node_iter == attribute_of_node.end() 
              || (*attribute_of_node_iter).second != *attr_iter))) {
            (*(distance_from_object_with_attribute_to_node[parentId]))
              [*attr_iter]->setNextNearestObjectWithAttributeValue(
                            nearest_object);
            (*(distance_from_object_with_attribute_to_node[parentId]))
              [*attr_iter]->setDistanceToNextNearestObjectWithAttributeValue(
                            attribute_to_parent_distance_through_node);
            have_updated_next_nearest_from_node_nearest = true;
          }
        }
        attribute_to_node_next_nearest_distance
        = (*(distance_from_object_with_attribute_to_node[nodeId]))
        [*attr_iter]->getDistanceToNextNearestObjectWithAttributeValue();
        if (!have_updated_next_nearest_from_node_nearest 
            && attribute_to_node_next_nearest_distance != -1.0) {
          const ObjectT next_nearest_object
          = (*(distance_from_object_with_attribute_to_node[nodeId]))
            [*attr_iter]->getNextNearestObjectWithAttributeValue();
          if (next_nearest_object != current_nearest_object_to_parent) {
            attribute_to_parent_next_nearest_distance_through_node
              = attribute_to_node_next_nearest_distance 
                + tree.getDistanceToFather(nodeId);
            if (current_attribute_to_parent_next_nearest_distance == -1.0 ||
              ((seek_longest
                ? (attribute_to_parent_next_nearest_distance_through_node
                  > current_attribute_to_parent_next_nearest_distance)
                : (attribute_to_parent_next_nearest_distance_through_node
                  < current_attribute_to_parent_next_nearest_distance)) ||
                attribute_to_parent_next_nearest_distance_through_node
                == current_attribute_to_parent_next_nearest_distance &&
                nearest_object < current_next_nearest_object_to_parent) &&
              (!suppress_inparalogs ||
              attribute_of_node_iter == attribute_of_node.end()
              || (*attribute_of_node_iter).second != *attr_iter)) {
              (*(distance_from_object_with_attribute_to_node[parentId]))
                [*attr_iter]->setNextNearestObjectWithAttributeValue(
                                next_nearest_object);
              (*(distance_from_object_with_attribute_to_node[parentId]))
              [*attr_iter]->setDistanceToNextNearestObjectWithAttributeValue(
                    attribute_to_parent_next_nearest_distance_through_node);
            }
          }
        }
      }
    }
  }
  if (seek_longest) {
    return;
  }
  int childId;
  double attribute_to_child_distance_through_node; 
  double current_attribute_to_child_distance, 
          attribute_to_child_next_nearest_distance_through_node,
          current_attribute_to_child_next_nearest_distance;
  for (int i = 0; i < breadth_first_visit_order.size(); i++) {
    nodeId = breadth_first_visit_order[i];
    const vector<int> &children = tree.getSonsId(breadth_first_visit_order[i]);
    for (int j = 0; j < children.size(); j++) {
      childId = children[j];
      for (attr_iter = attributes_with_multiple_objects.begin(); 
            attr_iter != attributes_with_multiple_objects.end(); 
              ++attr_iter) {
        have_updated_next_nearest_from_node_nearest = false;
        attribute_to_node_distance
          = (*(distance_from_object_with_attribute_to_node[nodeId]))
            [*attr_iter]->getDistanceToNearestObjectWithAttributeValue();
        if (attribute_to_node_distance != -1.0) {
          attribute_to_child_distance_through_node 
            = attribute_to_node_distance + tree.getDistanceToFather(childId);
          current_attribute_to_child_distance
            = (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->getDistanceToNearestObjectWithAttributeValue();
          const ObjectT nearest_object
            = (*(distance_from_object_with_attribute_to_node[nodeId]))
              [*attr_iter]->getNearestObjectWithAttributeValue();
          ObjectT current_nearest_object_to_child
            = (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->getNearestObjectWithAttributeValue();
          const ObjectT current_next_nearest_object_to_child
            = (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->getNextNearestObjectWithAttributeValue();
          current_attribute_to_child_next_nearest_distance
            = (*(distance_from_object_with_attribute_to_node[childId]))
            [*attr_iter]->getDistanceToNextNearestObjectWithAttributeValue();
          if (current_attribute_to_child_distance == -1.0 ||
                (seek_longest
                ? (attribute_to_child_distance_through_node
                    > current_attribute_to_child_distance)
                : (attribute_to_child_distance_through_node 
                    < current_attribute_to_child_distance)) ||
                attribute_to_child_distance_through_node
                  == current_attribute_to_child_distance &&
                (current_nearest_object_to_child == nullObjectValue ||
                  nearest_object < current_nearest_object_to_child)) {
            if (current_attribute_to_child_distance != -1.0 &&
                current_nearest_object_to_child != nullObjectValue &&
                current_nearest_object_to_child != nearest_object) {
              (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setNextNearestObjectWithAttributeValue(
                              current_nearest_object_to_child);
              (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setDistanceToNextNearestObjectWithAttributeValue(
                              current_attribute_to_child_distance);
              current_attribute_to_child_next_nearest_distance
                = current_attribute_to_child_distance;
            }
            (*(distance_from_object_with_attribute_to_node[childId]))
            [*attr_iter]->setNearestObjectWithAttributeValue(nearest_object);
            current_nearest_object_to_child = nearest_object;
            (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setDistanceToNearestObjectWithAttributeValue(
                attribute_to_child_distance_through_node);
            (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setNearestObjectIsDescendant(false);
          } else {
            if (nearest_object != current_nearest_object_to_child &&
                (current_attribute_to_child_next_nearest_distance == -1.0 ||
                ((seek_longest
                  ? (attribute_to_child_distance_through_node
                      < current_attribute_to_child_next_nearest_distance)
                  : (attribute_to_child_distance_through_node
                      < current_attribute_to_child_next_nearest_distance)) ||
                  attribute_to_child_distance_through_node
                  == current_attribute_to_child_next_nearest_distance &&
                  nearest_object < current_next_nearest_object_to_child))) {
              (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setNextNearestObjectWithAttributeValue(
                              nearest_object);
              (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setDistanceToNextNearestObjectWithAttributeValue(
                              attribute_to_child_distance_through_node);
              have_updated_next_nearest_from_node_nearest = true;
            }
          }
          attribute_to_node_next_nearest_distance
            = (*(distance_from_object_with_attribute_to_node[nodeId]))
            [*attr_iter]->getDistanceToNextNearestObjectWithAttributeValue();
          if (!have_updated_next_nearest_from_node_nearest && 
              attribute_to_node_next_nearest_distance != -1.0) {
            const ObjectT next_nearest_object
            = (*(distance_from_object_with_attribute_to_node[nodeId]))
              [*attr_iter]->getNextNearestObjectWithAttributeValue();
            if (next_nearest_object != current_nearest_object_to_child) {
              attribute_to_child_next_nearest_distance_through_node
                = attribute_to_node_next_nearest_distance 
                  + tree.getDistanceToFather(childId);
              if (current_attribute_to_child_next_nearest_distance == -1.0 ||
                ((seek_longest
                  ? (attribute_to_child_next_nearest_distance_through_node
                    > current_attribute_to_child_next_nearest_distance)
                  : (attribute_to_child_next_nearest_distance_through_node
                    < current_attribute_to_child_next_nearest_distance)) ||
                  attribute_to_child_next_nearest_distance_through_node
                  == current_attribute_to_child_next_nearest_distance &&
                  nearest_object < current_next_nearest_object_to_child)) {
                (*(distance_from_object_with_attribute_to_node[childId]))
                  [*attr_iter]->setNextNearestObjectWithAttributeValue(
                                  next_nearest_object);
                (*(distance_from_object_with_attribute_to_node[childId]))
              [*attr_iter]->setDistanceToNextNearestObjectWithAttributeValue(
                      attribute_to_child_next_nearest_distance_through_node);
              }
            }
          }
        }
      }
    }
  }
}
