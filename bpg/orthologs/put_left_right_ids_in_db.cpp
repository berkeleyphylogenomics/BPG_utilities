// put_left_right_ids_in_db.cpp
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

#include "put_left_right_ids_in_db.h"
#include <mysql++/query.h>
#include <mysql++/result.h>
#include <boost/regex.hpp>

using namespace std;
using namespace boost;

bool putLeftRightIdsInDB(const bpp::Tree &tree, 
                  const LeftRightIdsOfNodeIdMap &idMap,
                  const map<int, int> &nodeIdOfLeftIdMap,
                  const map<int, int> &levelOfNodeIdMap,
                  const DBIdentifierT tree_id,
                  const string &alignment_id,
                  mysqlpp::Connection &conn,
                  map<int, DBIdentifierT> &db_id_of_left_id,
                  map<int, DBIdentifierT> &sequence_header_of_leaf) {
  sequence_header_of_leaf.clear();
    mysqlpp::Query current_node_query = conn.query();
    current_node_query << "SELECT id, left_id, "
                       << "protein_sequence_id AS sequence_header_id "
                        << "FROM tree_node "
                        << "WHERE tree_id = " << tree_id;
    mysqlpp::ResUse current_node_res = current_node_query.use();
    int nodeId;
    DBIdentifierT protein_seq_id;
    bool alreadyInDatabase = false;
    map<int, int>::const_iterator node_id_of_left_id_iter;
    while (mysqlpp::Row row = current_node_res.fetch_row()) {
      alreadyInDatabase = true;
      db_id_of_left_id[row["left_id"]] = row["id"];
      protein_seq_id = int(row["sequence_header_id"]);
      if (0 != protein_seq_id) {
        node_id_of_left_id_iter = nodeIdOfLeftIdMap.find(row["left_id"]);
        nodeId = (*node_id_of_left_id_iter).second;
        sequence_header_of_leaf[nodeId] = protein_seq_id;
      }
    }
    if (alreadyInDatabase) {
      return false;
    }
    mysqlpp::Query node_query = conn.query();
    node_query  << "INSERT INTO tree_node "
                << "(tree_id, left_id, right_id, "
                << "parent_left_id, level) "
                << " VALUES "
                << "('" << tree_id << "', "
                << "%0q, %1q, %2q, %3q) "
                << "ON DUPLICATE KEY UPDATE "
                << "right_id = %1q, "
                << "parent_left_id = %2q, "
                << "level = %3q ";
    node_query.parse();
    vector<int> nodeIds = tree.getNodesId();
    vector<int>::const_iterator node_iter;
    LeftRightIdsOfNodeIdMap::const_iterator left_iter;
    map<int, int>::const_iterator level_iter;
    int root_id = tree.getRootId();
    int left_id, right_id, parent_id, parent_left_id, level;
    for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) {
      if (*node_iter == root_id) {
        parent_left_id = 0;
      } else {
        parent_id = tree.getFatherId(*node_iter);
        left_iter = idMap.find(parent_id);
        parent_left_id = ((*left_iter).second).first;
        level_iter = levelOfNodeIdMap.find(*node_iter);
        level = (*level_iter).second;
      }
      left_iter = idMap.find(*node_iter);
      left_id = ((*left_iter).second).first;
      right_id = ((*left_iter).second).second;
      mysqlpp::ResNSel status
        = node_query.execute(left_id, right_id, parent_left_id, level);
      if (conn.errnum()) {
        cerr << "Error: " << conn.error() << endl;
      }
    }
    mysqlpp::Query nodeDBIdQuery = conn.query();
    nodeDBIdQuery << "SELECT id, left_id FROM tree_node "
                << "WHERE tree_id = " << tree_id;
    mysqlpp::ResUse node_db_id_res = nodeDBIdQuery.use();
    while (mysqlpp::Row node_db_id_row = node_db_id_res.fetch_row()) {
      db_id_of_left_id[int(node_db_id_row["left_id"])] = node_db_id_row["id"];
    }
    vector<int> leafIds = tree.getLeavesId();
    map<string, vector<int> *> tree_node_ids_of_sequence_id;
    regex bpgseq_patt("bpgseq([0-9]+)");
    smatch seq_matches;
    vector<string> sequence_ids;
    string sequence_ids_clause = "(";
    for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
      if (regex_match(tree.getNodeName(*node_iter), seq_matches, bpgseq_patt)){
        if (tree_node_ids_of_sequence_id.find(seq_matches[1])
            == tree_node_ids_of_sequence_id.end()) {
          tree_node_ids_of_sequence_id[seq_matches[1]] = new vector<int>();
        }
        tree_node_ids_of_sequence_id[seq_matches[1]]->push_back(*node_iter);
        sequence_ids.push_back(seq_matches[1]);
        sequence_ids_clause += seq_matches[1];
        sequence_ids_clause += ",";
      }
      sequence_header_of_leaf[*node_iter] = 0;
    }
    // Since we have an extra comma at the end, we append a nonexistent sequence
    // id of 0 to make the SQL clause well-formed
    sequence_ids_clause += "0) ";
    mysqlpp::Query sequence_query = conn.query();
    sequence_query << "SELECT sequence_id, "
                  << "protein_sequence_id AS sequence_header_id "
                  << "FROM alignment_protein_sequence, protein_sequence "
                  << "WHERE alignment_id = " << alignment_id
                  << " AND protein_sequence_id = protein_sequence.id "
                  << "AND sequence_id IN " << sequence_ids_clause
                  << "ORDER BY sequence_id";
    mysqlpp::ResUse sequence_res = sequence_query.use();
    int i = 0, leafId;
    string seq_id = "";
    protein_seq_id = 0;
    while (mysqlpp::Row row = sequence_res.fetch_row()) {
      if (row["sequence_id"].get_string() == seq_id) {
        if (i < tree_node_ids_of_sequence_id[seq_id]->size() - 1) {
          ++i;
        }
      } else {
        if (seq_id != "" && protein_seq_id != 0) {
          ++i;
          while (i < tree_node_ids_of_sequence_id[seq_id]->size()) {
            leafId = (*tree_node_ids_of_sequence_id[seq_id])[i];
            sequence_header_of_leaf[leafId] = protein_seq_id;
            ++i;
          }
        }
        i = 0;
      }
      seq_id = row["sequence_id"].get_string();
      protein_seq_id = int(row["sequence_header_id"]);
      leafId = (*tree_node_ids_of_sequence_id[seq_id])[i];
      sequence_header_of_leaf[leafId] = protein_seq_id;
    }
    mysqlpp::Query leaf_query = conn.query();
    leaf_query  << "UPDATE tree_node "
                << "SET protein_sequence_id = %0q "
                << "WHERE tree_id = " << tree_id
                << " AND left_id = %1q";
    leaf_query.parse();
    map<int, DBIdentifierT>::const_iterator sequence_header_of_leaf_iter;
    for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) {
      sequence_header_of_leaf_iter = sequence_header_of_leaf.find(*node_iter);
      if (sequence_header_of_leaf_iter == sequence_header_of_leaf.end()) {
        cerr << "No sequence header assigned to node " << *node_iter << endl;
      } else {
        left_iter = idMap.find(*node_iter);
        left_id = ((*left_iter).second).first;
        mysqlpp::ResNSel status
          = leaf_query.execute(sequence_header_of_leaf[*node_iter], left_id);
      }
    }
    map<string, vector<int> *>::iterator node_of_seq_iter;
    for(node_of_seq_iter = tree_node_ids_of_sequence_id.begin();
        node_of_seq_iter != tree_node_ids_of_sequence_id.end();
        ++node_of_seq_iter) {
      delete (*node_of_seq_iter).second;
    }
  return true;
}
