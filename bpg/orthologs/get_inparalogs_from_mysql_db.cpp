// get_inparalogs_from_mysql_db.cpp
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

#include "get_inparalogs_from_db.h"
#include <mysql++/connection.h>
#include <mysql++/query.h>
#include <mysql++/result.h>

//template<class AttributeT>
void getInparalogsFromDB(const bpp::Tree &tree,
                        const map<int, int> &nodeIdOfLeftIdMap,
                        const int tree_id,
                        const string &attribute_name,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
                        map<int, AttributeT > &attribute_of_node) {
  attribute_of_node.clear();
    mysqlpp::Query inparalog_query = conn.query();
    inparalog_query << "SELECT left_id, " << attribute_name 
                << " FROM tree_node, " 
                << "tree_node_" << attribute_name
                << " WHERE tree_id = " << tree_id
                << " AND tree_node.id = tree_node_id";
    mysqlpp::ResUse inparalog_res = inparalog_query.use();
    int nodeId;
    while (mysqlpp::Row row = inparalog_res.fetch_row()) {
      nodeId = (*(nodeIdOfLeftIdMap.find(row["left_id"]))).second;
      attribute_of_node[nodeId] = int(row[attribute_name.c_str()]);
    }
}