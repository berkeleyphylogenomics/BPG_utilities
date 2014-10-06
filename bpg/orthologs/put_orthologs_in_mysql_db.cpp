// put_orthologs_in_mysql_db.cpp
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

#include "put_orthologs_in_db.h"
#include <mysql++/query.h>
#include <mysql++/result.h>
#include <iostream>

//template<class ObjectT, ObjectT nullObjectValue,
//        class AttributeT, AttributeT nullAttributeValue>
bool putOrthologsInDB(const bpp::Tree &tree,
            const LeftRightIdsOfNodeIdMap &idMap,
            const map<int, DBIdentifierT> &db_id_of_left_id,
            const string &attribute_name,
            const string &object_name,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
            const map<AttributeT , ObjectT > 
              unique_object_with_attribute,
            const map<int, map<AttributeT , TreeDistanceInfo<ObjectT, 
                                                        nullObjectValue> *> *>
              distance_from_object_with_attribute_to_node,
            const map<int, set<AttributeT > *> 
              &attributes_with_multiple_objects_for_which_node_is_maximal) {
    LeftRightIdsOfNodeIdMap::const_iterator left_iter;
    map<int, DBIdentifierT>::const_iterator db_left_iter;
    int rootId = tree.getRootId();
    left_iter = idMap.find(rootId);
    int leftIdOfRootId = ((*left_iter).second).first;
    db_left_iter = db_id_of_left_id.find(leftIdOfRootId);
    int dbIdOfRootId = (*db_left_iter).second;
    int dbIdToTest;
//    typename 
    map<int, set<AttributeT > *>::const_iterator node_attrs_iter;
    if (unique_object_with_attribute.size() > 0) {
      dbIdToTest = dbIdOfRootId;
    } else if (attributes_with_multiple_objects_for_which_node_is_maximal.size()
                > 0) {
      node_attrs_iter 
        = attributes_with_multiple_objects_for_which_node_is_maximal.begin();
      int nodeId = (*node_attrs_iter).first;
      left_iter = idMap.find(nodeId);
      int leftIdOfNodeId = ((*left_iter).second).first;
      db_left_iter = db_id_of_left_id.find(leftIdOfNodeId);
      dbIdToTest = (*db_left_iter).second;
    } else {
      // Nothing to do
      return false;
    }
    mysqlpp::Query current_ortholog_query = conn.query();
    current_ortholog_query 
        << "SELECT tree_node_id "
        << "FROM "
        << "maximal_tree_node_nearest_" << object_name
        << "_with_" << attribute_name
        << " WHERE tree_node_id = " << dbIdToTest;
    mysqlpp::ResUse current_ortholog_res = current_ortholog_query.use();
    bool alreadyInDatabase = false;
    while (mysqlpp::Row row = current_ortholog_res.fetch_row()) {
      alreadyInDatabase = true;
    }
    if (alreadyInDatabase) {
      return false;
    }
    mysqlpp::Query ortholog_query = conn.query();
    ortholog_query << "INSERT INTO "
                << "maximal_tree_node_nearest_" << object_name
                << "_with_" << attribute_name
                << " (tree_node_id, "
                << attribute_name << ", "
                << object_name << ", "
                << "distance, "
                << "next_nearest_" << object_name << ", "
                << "next_nearest_distance, is_descendant, is_unique) "
                << " VALUES "
                << "(%0q, %1q, %2q, %3q, %4q, %5q, %6q, %7q) "
                << "ON DUPLICATE KEY UPDATE "
                << object_name << " = '%2q', "
                << "distance = %3q, "
                << "next_nearest_" << object_name << " = '%4q', "
                << "next_nearest_distance = %5q, "
                << "is_descendant = %6q, is_unique = %7q ";
    ortholog_query.parse();
//    typename 
    map<AttributeT , ObjectT >::const_iterator 
        obj_attr_iter;
    for (obj_attr_iter = unique_object_with_attribute.begin();
          obj_attr_iter != unique_object_with_attribute.end();
          ++obj_attr_iter) {
      const AttributeT attr = (*obj_attr_iter).first;
      const ObjectT obj = (*obj_attr_iter).second;
      mysqlpp::ResNSel status
        = ortholog_query.execute(dbIdOfRootId, attr, obj, -1.0, 
                                nullObjectValue, -1.0, true, true);
      if (conn.errnum()) {
        cerr << "Error: " << conn.error() << endl;
      }
    }
//    typename 
    set<AttributeT >::const_iterator attr_iter;
//    typename 
    map<int, map<AttributeT , TreeDistanceInfo<ObjectT,
                    nullObjectValue> *> *>::const_iterator dist_attr_node_iter;
//    typename 
    map<AttributeT, TreeDistanceInfo<ObjectT, nullObjectValue> *>::const_iterator 
          attr_node_iter;
    int nodeId, leftId, dbId;
    double distance, next_nearest_distance;
    bool isDescendant;
    for (node_attrs_iter = 
          attributes_with_multiple_objects_for_which_node_is_maximal.begin();
        node_attrs_iter !=
          attributes_with_multiple_objects_for_which_node_is_maximal.end();
        ++node_attrs_iter) {
      nodeId = (*node_attrs_iter).first;
      left_iter = idMap.find(nodeId);
      leftId = ((*left_iter).second).first;
      db_left_iter = db_id_of_left_id.find(leftId);
      dbId = (*db_left_iter).second;
      for (attr_iter = (*node_attrs_iter).second->begin();
          attr_iter != (*node_attrs_iter).second->end(); ++attr_iter) {
        dist_attr_node_iter 
          = distance_from_object_with_attribute_to_node.find(nodeId);
        attr_node_iter = ((*dist_attr_node_iter).second)->find(*attr_iter);
        const ObjectT obj
          = ((*attr_node_iter).second)->getNearestObjectWithAttributeValue();
        const ObjectT next_nearest_obj
        = ((*attr_node_iter).second)->getNextNearestObjectWithAttributeValue();
        distance = ((*attr_node_iter).second)
                    ->getDistanceToNearestObjectWithAttributeValue();
        next_nearest_distance = ((*attr_node_iter).second)
                    ->getDistanceToNextNearestObjectWithAttributeValue();
        isDescendant = ((*attr_node_iter).second)->isNearestObjectDescendant();
        mysqlpp::ResNSel status
          = ortholog_query.execute(dbId, *attr_iter, obj, distance, 
                  next_nearest_obj, next_nearest_distance, isDescendant, false);
      }
    }
  return true;
}
