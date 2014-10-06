// put_superorthologous_nodes_in_postgresql_db.cpp
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

#include "put_superorthologous_nodes_in_db.h"
#include <pqxx/transactor>
#include <iostream>

void getContainedLeaves(const bpp::Tree &tree, int nodeId, 
                        vector<int> &leaves) {
  if (tree.isLeaf(nodeId)) {
    leaves.push_back(nodeId);
  } else {
    const vector<int> &children = tree.getSonsId(nodeId);
    for (vector<int>::const_iterator child_iter = children.begin();
        child_iter != children.end(); ++child_iter) {
      getContainedLeaves(tree, *child_iter, leaves);
    }
  }
}
        
void putSuperorthologousNodesInDB(const bpp::Tree &tree,
            const LeftRightIdsOfNodeIdMap &idMap,
            const map<int, DBIdentifierT> &db_id_of_left_id,
            const int tree_id,
            const string &attribute_name,
            const string &object_name,
            pqxx::connection &conn,
            const set<int> &super_orthologous_nodes) {
    LeftRightIdsOfNodeIdMap::const_iterator left_iter;
    int left_id, right_id;
    DBIdentifierT db_id;
    map<int, DBIdentifierT>::const_iterator db_left_iter;
    set<int>::const_iterator node_iter;
    vector<int> leaves;
    pqxx::work T(conn, "putSuperorthologousNodesInDBTransaction");
    ostringstream reset_query;
    reset_query << "UPDATE tree_node "
            << "SET is_superorthologous = FALSE "
            << "WHERE "
            << "tree_id = " << tree_id
            << " AND is_superorthologous = TRUE";
    T.exec(reset_query.str());
    for (node_iter = super_orthologous_nodes.begin();
          node_iter != super_orthologous_nodes.end(); ++node_iter) {
      left_iter = idMap.find(*node_iter);
      left_id = ((*left_iter).second).first;
      right_id = ((*left_iter).second).second;
      db_left_iter = db_id_of_left_id.find(left_id);
      db_id = (*db_left_iter).second;
      ostringstream node_query;
      node_query << "UPDATE tree_node "
              << "SET is_superorthologous = TRUE "
              << "WHERE "
              << "tree_id = " << tree_id
              << " AND left_id = " << left_id;
      T.exec(node_query.str());
      leaves.clear();
      getContainedLeaves(tree, *node_iter, leaves);
      ostringstream leftid_clause;
      leftid_clause << "(";
      for (vector<int>::const_iterator leaf_iter = leaves.begin();
            leaf_iter != leaves.end(); ++leaf_iter) {
        left_iter = idMap.find(*leaf_iter);
        leftid_clause << left_iter->second.first << ",";
      }
      leftid_clause << "0)";
      ostringstream superortholog_query;
      superortholog_query << "UPDATE tree_node "
              << "SET superorthologous_node_id = " << db_id << ", "
              << "superorthologous_node_left_id = " << left_id
              << " WHERE tree_id = " << tree_id
              << " AND left_id IN " << leftid_clause.str();
      T.exec(superortholog_query.str());
    }
    T.commit();
}
