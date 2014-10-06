// get_orthologs_from_mysql_db.cpp
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

#include "get_orthologs_from_db.h"
#include <mysql++/query.h>
#include <mysql++/result.h>

//template<class ObjectT, ObjectT nullObjectValue,
//        class AttributeT, AttributeT nullAttributeValue>
void getOrthologsFromDB(const bpp::Tree &tree,
            const map<int, int> &nodeIdOfLeftIdMap,
            const int tree_id,
            const string &attribute_name,
            const string &object_name,
            mysqlpp::Connection &conn,
            map<AttributeT , ObjectT > 
                &unique_object_with_attribute,
            map<int, map<AttributeT , TreeDistanceInfo<ObjectT, 
                                                        nullObjectValue> *> *>
                &distance_from_object_with_attribute_to_node,
            map<int, set<AttributeT > *> 
                &attributes_with_multiple_objects_for_which_node_is_maximal) {
  unique_object_with_attribute.clear();
  distance_from_object_with_attribute_to_node.clear();
  attributes_with_multiple_objects_for_which_node_is_maximal.clear();
    mysqlpp::Query ortholog_query = conn.query();
    ortholog_query << "SELECT left_id, "
                    << attribute_name << ", "
                    << "maximal_tree_node_nearest_" << object_name
                    << "_with_" << attribute_name
                    << "." << object_name << ", "
                    << "distance, "
                    << "next_nearest_" << object_name << ", "
                    << "next_nearest_distance, is_descendant, "
                    << "is_unique "
                    << "FROM tree_node, "
                    << "maximal_tree_node_nearest_" << object_name
                    << "_with_" << attribute_name
                    << " WHERE tree_id = " << tree_id
                    << " AND tree_node.id = tree_node_id";
    mysqlpp::ResUse ortholog_res = ortholog_query.use();
    if (conn.errnum()) {
      cerr << "Error: " << conn.error() << endl;
    }
    int nodeId;
    bool isUnique, isDescendant;
    double distance, next_nearest_distance;
    while (mysqlpp::Row row = ortholog_res.fetch_row()) {
      nodeId = (*(nodeIdOfLeftIdMap.find(row["left_id"]))).second;
      AttributeT attr = int(row[attribute_name.c_str()]);
      ObjectT obj = int(row[object_name.c_str()]);
      distance = double(row["distance"]);
      next_nearest_distance = double(row["next_nearest_distance"]);
      isUnique = int(row["is_unique"]);
      isDescendant = int(row["is_descendant"]);
      if (isUnique) {
        unique_object_with_attribute[attr] = obj;
      } else {
        if (attributes_with_multiple_objects_for_which_node_is_maximal.find(
                nodeId) ==
            attributes_with_multiple_objects_for_which_node_is_maximal.end()) {
          attributes_with_multiple_objects_for_which_node_is_maximal[nodeId]
            = new set<AttributeT>;
        }
        attributes_with_multiple_objects_for_which_node_is_maximal[
                                                        nodeId]->insert(attr);
        if (distance_from_object_with_attribute_to_node.find(nodeId) 
            == distance_from_object_with_attribute_to_node.end()) {
          distance_from_object_with_attribute_to_node[nodeId]
            = new map<AttributeT , TreeDistanceInfo<ObjectT, 
                                                nullObjectValue> *>;
        }
        if (distance_from_object_with_attribute_to_node[nodeId]->find(attr)
            == distance_from_object_with_attribute_to_node[nodeId]->end()) {
          (*(distance_from_object_with_attribute_to_node[nodeId]))[attr]
            = new TreeDistanceInfo<ObjectT, nullObjectValue>();
        }
        (*(distance_from_object_with_attribute_to_node[nodeId]))
          [attr]->setNearestObjectWithAttributeValue(obj);
        (*(distance_from_object_with_attribute_to_node[nodeId]))
          [attr]->setDistanceToNearestObjectWithAttributeValue(distance);
        (*(distance_from_object_with_attribute_to_node[nodeId]))
          [attr]->setDistanceToNextNearestObjectWithAttributeValue(
                                                    next_nearest_distance);
        (*(distance_from_object_with_attribute_to_node[nodeId]))
          [attr]->setNearestObjectIsDescendant(isDescendant);
          
      }
    }
}
