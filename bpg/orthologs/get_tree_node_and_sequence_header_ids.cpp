// get_tree_node_and_sequence_header_ids.cpp
// Author: Ruchira S. Datta
// Copyright (c) 2010, Regents of the University of California
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

#include "get_tree_node_and_sequence_header_ids.h"
#include <pqxx/transactor>
#include <pqxx/result>
#include <boost/regex.hpp>

using namespace std;
using namespace boost;

bool getTreeNodeAndSequenceHeaderIds(const bpp::Tree &tree, 
                  const LeftRightIdsOfNodeIdMap &idMap,
                  const map<int, int> &nodeIdOfLeftIdMap,
                  const map<int, int> &levelOfNodeIdMap,
                  const DBIdentifierT tree_id,
                  pqxx::connection &conn,
                  map<int, DBIdentifierT> &db_id_of_left_id,
                  map<int, DBIdentifierT> &sequence_header_of_leaf,
                  map<DBIdentifierT, int> &leaf_of_sequence_header) {
  db_id_of_left_id.clear();
  sequence_header_of_leaf.clear();
  ostringstream tree_node_query;
  tree_node_query << "SELECT id, left_id, sequence_header_id "
                  << "FROM tree_node "
                  << "WHERE tree_id = " << tree_id;
  pqxx::result R;
  pqxx::work T(conn, "getTreeNodeAndSequenceHeaderIdsTransaction");
  R = T.exec(tree_node_query.str());
  DBIdentifierT db_id;
  int left_id;
  int nodeId;
  DBIdentifierT sequence_header_id;
  map<int, int>::const_iterator node_id_of_left_id_iter;
  for (pqxx::result::const_iterator i = R.begin(); i != R.end(); ++i) {
    (*i)[0].to(db_id);
    (*i)[1].to(left_id);
    db_id_of_left_id[left_id] = db_id;
    if (!(*i)[2].is_null()) {
      (*i)[2].to(sequence_header_id);
        node_id_of_left_id_iter = nodeIdOfLeftIdMap.find(left_id);
        nodeId = (*node_id_of_left_id_iter).second;
        sequence_header_of_leaf[nodeId] = sequence_header_id;
        leaf_of_sequence_header[sequence_header_id] = nodeId;
    }
  }
  return true;
}
