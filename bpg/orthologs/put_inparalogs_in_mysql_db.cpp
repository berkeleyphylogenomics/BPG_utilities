// put_inparalogs_in_mysql_db.cpp
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


#include "put_inparalogs_in_db.h"
#include <mysql++/connection.h>
#include <mysql++/query.h>
#include <mysql++/result.h>

//template<class AttributeT>
bool putInparalogsInDB(const bpp::Tree &tree,
                        const LeftRightIdsOfNodeIdMap &idMap,
                        const map<int, DBIdentifierT> &db_id_of_left_id,
                        const string &attribute_name,
                        mysqlpp::Connection &conn,
                        const map<int, AttributeT > &attribute_of_node) {
    LeftRightIdsOfNodeIdMap::const_iterator left_iter;
    ostringstream tree_node_ids_clause;
    tree_node_ids_clause << "(";
    map<int, AttributeT >::const_iterator attr_of_node_iter;
    map<int, DBIdentifierT>::const_iterator db_left_iter;
    int nodeId, leftId, dbId;
    for (attr_of_node_iter = attribute_of_node.begin();
          attr_of_node_iter != attribute_of_node.end();
          ++attr_of_node_iter) {
      nodeId = (*attr_of_node_iter).first;
      left_iter = idMap.find(nodeId);
      leftId = ((*left_iter).second).first;
      db_left_iter = db_id_of_left_id.find(leftId);
      dbId = (*db_left_iter).second;
      tree_node_ids_clause << dbId << ",";
    }
    tree_node_ids_clause << "0)";
    mysqlpp::Query current_attribute_of_node_query = conn.query();
    current_attribute_of_node_query 
        << "SELECT tree_node_id "
        << "FROM " << "tree_node_" << attribute_name
        << " WHERE tree_node_id "
        << " IN " << tree_node_ids_clause.str()
        << " LIMIT 1";
    mysqlpp::ResUse current_res = current_attribute_of_node_query.use();
    bool alreadyInDatabase = false;
    while (mysqlpp::Row current_row = current_res.fetch_row()) {
      alreadyInDatabase = true;
    }
    if (alreadyInDatabase) {
      return false;
    }
    mysqlpp::Query attribute_of_node_query = conn.query();
    attribute_of_node_query << "INSERT INTO " << "tree_node_" << attribute_name
                            << " (tree_node_id, "
                            << attribute_name << ", "
                            << "is_maximal) "
                            << " VALUES "
                            << "(%0q, %1q, %2q) "
                            << "ON DUPLICATE KEY UPDATE "
                            << "is_maximal = %2q ";
    attribute_of_node_query.parse();
//      typename 
    int rootId = tree.getRootId();
    bool isMaximal;
    for (attr_of_node_iter = attribute_of_node.begin();
        attr_of_node_iter != attribute_of_node.end();
        ++attr_of_node_iter) {
      nodeId = (*attr_of_node_iter).first;
      AttributeT attr = (*attr_of_node_iter).second;
      isMaximal = true;
      if (nodeId != rootId && 
          attribute_of_node.find(tree.getFatherId(nodeId)) 
            != attribute_of_node.end()) {
        isMaximal = false;
      }
      left_iter = idMap.find(nodeId);
      leftId = ((*left_iter).second).first;
      db_left_iter = db_id_of_left_id.find(leftId);
      dbId = (*db_left_iter).second;
      mysqlpp::ResNSel status
        = attribute_of_node_query.execute(dbId, attr, isMaximal);
      if (conn.errnum()) {
        cerr << "Error: " << conn.error() << endl;
      }
    }
  return true;
}