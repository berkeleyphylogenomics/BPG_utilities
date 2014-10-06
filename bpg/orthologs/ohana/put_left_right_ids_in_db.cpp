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
#include <boost/regex.hpp>
#include <sstream>
#include <iostream>

using namespace std;
using namespace boost;

bool putLeftRightIdsInDB(const bpp::Tree &tree, 
                  const LeftRightIdsOfNodeIdMap &idMap,
                  const map<int, int> &nodeIdOfLeftIdMap,
                  const map<int, int> &levelOfNodeIdMap,
                  const int tree_id,
                  const string &alignment_id,
                  map<int, int> &db_id_of_left_id,
                  map<int, int> &protein_sequence_of_leaf,
                  pqxx::transaction<>& xaction) 
{
    protein_sequence_of_leaf.clear();

    int nodeId;
    int protein_seq_id;
    map<int, int>::const_iterator node_id_of_left_id_iter;
    vector<int> nodeIds = tree.getNodesId();
    vector<int>::const_iterator node_iter;
    LeftRightIdsOfNodeIdMap::const_iterator left_iter;
    map<int, int>::const_iterator level_iter;
    int root_id = tree.getRootId();
    int left_id, right_id, parent_id, parent_left_id, level;

    for (node_iter = nodeIds.begin(); node_iter != nodeIds.end(); ++node_iter) 
    {
        if (*node_iter == root_id)
            parent_left_id = 0;
         
        else 
        {
            parent_id = tree.getFatherId(*node_iter);
            left_iter = idMap.find(parent_id);
            parent_left_id = ((*left_iter).second).first;
            level_iter = levelOfNodeIdMap.find(*node_iter);
            level = (*level_iter).second;
        }

        left_iter = idMap.find(*node_iter);
        left_id = ((*left_iter).second).first;
        right_id = ((*left_iter).second).second;

        ostringstream query;

        query << "INSERT INTO tree_node (tree_id, left_id, right_id, parent_left_id, level) VALUES ("
              << tree_id << ", " << left_id << ", " << right_id << ", " << parent_left_id
              << ", " << level << ")";

        xaction.exec(query.str());
    }

    ostringstream query;

    query << "SELECT id, left_id FROM tree_node WHERE tree_id=" << tree_id;
    pqxx::result res = xaction.exec(query.str());

    for(int i=0; i<res.size(); i++)
        db_id_of_left_id[atoi(res[i]["left_id"].c_str())] = atoi(res[i]["id"].c_str());

    vector<int> leafIds = tree.getLeavesId();
    map<string, vector<int> *> tree_node_ids_of_sequence_id;
    regex bpgseq_patt("seqh([0-9]+)");
    smatch seq_matches;
    vector<string> sequence_ids;
    string sequence_ids_clause = "(";

    for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) 
    {
        if (regex_match(tree.getNodeName(*node_iter), seq_matches, bpgseq_patt))
        {
            if (tree_node_ids_of_sequence_id.find(seq_matches[1]) == tree_node_ids_of_sequence_id.end()) 
            {
                tree_node_ids_of_sequence_id[seq_matches[1]] = new vector<int>();
            }

            tree_node_ids_of_sequence_id[seq_matches[1]]->push_back(*node_iter);
            sequence_ids.push_back(seq_matches[1]);
            sequence_ids_clause += seq_matches[1];
            sequence_ids_clause += ",";
        }

        protein_sequence_of_leaf[*node_iter] = 0;
    }

    // Since we have an extra comma at the end, we append a nonexistent sequence
    // id of 0 to make the SQL clause well-formed
    sequence_ids_clause += "0) ";

    query.str("");
    query << "SELECT sequence_header_id FROM family_sequence_header WHERE family_id="
          << alignment_id << " AND sequence_header_id IN " << sequence_ids_clause
          << " ORDER BY sequence_header_id";

    res = xaction.exec(query.str());
    
    int j = 0, leafId;
    string seq_id = "";
    protein_seq_id = 0;

    for(int i=0; i<res.size(); i++)
    {
        if(seq_id.compare(res[i]["sequence_header_id"].c_str()) == 0)
        {
            if(j < tree_node_ids_of_sequence_id[seq_id]->size() - 1)
                j++;
        }
    
        else
        {
            if(seq_id.compare("") != 0 && protein_seq_id != 0)
            {
                j++;

                while(j < tree_node_ids_of_sequence_id[seq_id]->size())
                {
                    leafId = (*tree_node_ids_of_sequence_id[seq_id])[j];
                    protein_sequence_of_leaf[leafId] = protein_seq_id;
                    j++;
                }
            }

            j = 0;
        }

        seq_id = res[i]["sequence_header_id"].c_str();
        protein_seq_id = atoi(seq_id.c_str());
        leafId = (*tree_node_ids_of_sequence_id[seq_id])[j];
        protein_sequence_of_leaf[leafId] = protein_seq_id;
    }
/*

    mysqlpp::Query leaf_query = conn.query();
    leaf_query  << "UPDATE tree_node "
                << "SET protein_sequence_id = %0q "
                << "WHERE tree_id = " << tree_id
                << " AND left_id = %1q";
    leaf_query.parse();
*/

    map<int, int>::const_iterator protein_sequence_of_leaf_iter;
    
    for (node_iter = leafIds.begin(); node_iter != leafIds.end(); ++node_iter) 
    {
        protein_sequence_of_leaf_iter = protein_sequence_of_leaf.find(*node_iter);
      
        if (protein_sequence_of_leaf_iter == protein_sequence_of_leaf.end()) 
        {
            cerr << "No protein sequence assigned to node " << *node_iter << endl;
        } 
        else 
        {
            left_iter = idMap.find(*node_iter);
            left_id = ((*left_iter).second).first;

            query.str("");
            query << "UPDATE tree_node SET sequence_id=" << protein_sequence_of_leaf[*node_iter]
                  << " WHERE tree_id=" << tree_id << " AND left_id=" << left_id;
            xaction.exec(query.str());        
        }
    }
    
    map<string, vector<int> *>::iterator node_of_seq_iter;
        
    for(node_of_seq_iter = tree_node_ids_of_sequence_id.begin(); node_of_seq_iter != tree_node_ids_of_sequence_id.end(); ++node_of_seq_iter) 
        delete (*node_of_seq_iter).second;

    return true;
}

