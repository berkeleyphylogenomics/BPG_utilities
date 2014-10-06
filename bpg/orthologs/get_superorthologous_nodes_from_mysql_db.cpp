// get_superorthologous_nodes_from_mysql_db.h
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

#include "get_superorthologous_nodes_from_db.h"
#include <mysql++/connection.h>
#include <mysql++/query.h>
#include <mysql++/result.h>

void getSuperorthologousNodesFromDB(const bpp::Tree &tree,
            const map<int, int> &nodeIdOfLeftIdMap,
            const int tree_id,
            const string &attribute_name,
            const string &object_name,
            mysqlpp::Connection &conn,
            set<int> &super_orthologous_nodes) {
  super_orthologous_nodes.clear();
    mysqlpp::Query superorthologous_query = conn.query();
    superorthologous_query << "SELECT left_id "
                          << "FROM tree_node "
                          << " WHERE tree_id = " << tree_id
                          << " AND is_superorthologous = TRUE";
    mysqlpp::ResUse superortholog_res = superorthologous_query.use();
    int nodeId;
    while (mysqlpp::Row row = superortholog_res.fetch_row()) {
      nodeId = (*(nodeIdOfLeftIdMap.find(row["left_id"]))).second;
      super_orthologous_nodes.insert(nodeId);
    }
}
