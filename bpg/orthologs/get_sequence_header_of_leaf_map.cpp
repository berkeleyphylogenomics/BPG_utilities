// get_sequence_header_of_leaf_map.cpp
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

#include "get_sequence_header_of_leaf_map.h"
#ifdef USING_POSTGRES
#include <pqxx/transactor>
#include <pqxx/result>
#include <iostream>
#else
#include <mysql++/connection.h>
#include <mysql++/query.h>
#include <mysql++/result.h>
#endif

void getSequenceHeaderOfLeafMap(const int tree_id,
                    const map<int, int> &nodeIdOfLeftIdMap,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
                    map<int, DBIdentifierT> 
                        &sequence_header_of_leaf) {
  sequence_header_of_leaf.clear();
#ifdef USING_POSTGRES
	pqxx::work T(conn, "getSequenceHeaderOfLeafMapTransaction");
	ostringstream query;
  query << "SELECT left_id, "
        << "sequence_header_id "
        << "FROM tree_node "
        << "WHERE tree_id = "
        << tree_id
        << " AND sequence_header_id IS NOT NULL "
        << "AND sequence_header_id <> 0 "
        << "ORDER BY left_id";
	pqxx::result R = T.exec(query.str());
	int left_id, leafId;
	DBIdentifierT sequence_header_id;
	for (pqxx::result::const_iterator i = R.begin(); i != R.end(); ++i) {
		(*i)[0].to(left_id);
		(*i)[1].to(sequence_header_id);
    leafId = (*(nodeIdOfLeftIdMap.find(left_id))).second;
		sequence_header_of_leaf[leafId] = sequence_header_id;
	}
#else
    mysqlpp::Query leaf_query = conn.query();
    leaf_query  << "SELECT left_id, "
                << "protein_sequence_id "
                << "FROM tree_node "
                << "WHERE tree_id = " << tree_id
                << " AND protein_sequence_id IS NOT NULL "
                << " AND protein_sequence_id <> 0 "
                << " ORDER BY left_id";
    mysqlpp::ResUse leaf_res = leaf_query.use();
    int leafId;
    while (mysqlpp::Row row = leaf_res.fetch_row()) {
      leafId = (*(nodeIdOfLeftIdMap.find(row["left_id"]))).second;
      sequence_header_of_leaf[leafId] = int(row["protein_sequence_id"]);
    }
#endif
}
