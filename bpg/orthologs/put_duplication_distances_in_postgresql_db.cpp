// put_duplication_distances_in_postgresql_db.cpp
// Author: Ruchira S. Datta
// Copyright (c) 2011, Regents of the University of California
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

#include "put_duplication_distances_in_db.h"
#include <pqxx/transactor>
#include <pqxx/result>
#include <iostream>
#include <assert.h>
using namespace std;

//template<class ObjectT, ObjectT nullObjectValue>
bool putDuplicationDistancesInDB(const bpp::Tree &tree,
    const LeftRightIdsOfNodeIdMap &idMap,
    const int tree_id,
    const string object_name,
    pqxx::connection &conn,
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *> &pair_of_nearest_objects_to_maximal_node,
    const map<int,  double> &greatest_distance_of_maximal_descendant) {
  pqxx::result R;
  pqxx::work T(conn, "putDuplicationDistancesInDBTransaction");
  map<int, double>::const_iterator max_desc_dist_iter;
  map<int, ObjObjInfo *>::const_iterator int_objobj_iter;
  double duplication_distance;
  int leftId;
  for (max_desc_dist_iter = greatest_distance_of_maximal_descendant.begin();
      max_desc_dist_iter != greatest_distance_of_maximal_descendant.end();
      ++max_desc_dist_iter) {
    leftId = idMap.find(max_desc_dist_iter->first)->second.first;
    ostringstream tree_node_query;
    int_objobj_iter 
      = pair_of_nearest_objects_to_maximal_node.find(max_desc_dist_iter->first);
    if (int_objobj_iter == pair_of_nearest_objects_to_maximal_node.end()) {
      // This had better be the root
      assert(leftId == 1);
      duplication_distance = 5000.0;
    } else {
      duplication_distance = getMaximalNodeDistance(max_desc_dist_iter->first,
                            duplication_node_of_objects,
                            pair_of_objects_yielding_duplication_node_distance,
                            pair_of_nearest_objects_to_maximal_node);

    }
    tree_node_query << "UPDATE tree_node "
        << "SET greatest_duplication_distance_of_maximal_descendant = "
        << max_desc_dist_iter->second
        << ", duplication_distance = " << duplication_distance
        << " WHERE "
        << "tree_id = " << tree_id
        << " AND left_id = " << leftId;
    T.exec(tree_node_query.str());
  }

  // We only insert the supporting evidence if it hasn't already been inserted
  bool alreadyInDatabase = false;
  ostringstream current_dup_evidence_query;
  current_dup_evidence_query << "SELECT left_id "
    << " FROM duplication_distance_of_duplication_node "
    << " WHERE tree_id = " << tree_id
    << " LIMIT 1";
  R = T.exec(current_dup_evidence_query.str());
  if (R.size() > 0) {
    alreadyInDatabase = true;
  } else {
    DupInfo *dupInfoPtr;
    ObjectT lesserObj, greaterObj;
    int otherLeftId;
    for (int_objobj_iter 
        = pair_of_objects_yielding_duplication_node_distance.begin();
        int_objobj_iter
        != pair_of_objects_yielding_duplication_node_distance.end();
        ++int_objobj_iter) {
      leftId = idMap.find(int_objobj_iter->first)->second.first;
      lesserObj = int_objobj_iter->second->getLesserObj();
      greaterObj = int_objobj_iter->second->getGreaterObj();
      dupInfoPtr = duplication_node_of_objects.find(lesserObj)->second->find(
                    greaterObj)->second;
      otherLeftId 
        = idMap.find(dupInfoPtr->getDuplicationNodeId())->second.first;
      ostringstream dup_node_query;
      duplication_distance = getDuplicationNodeDistance(int_objobj_iter->first,
                          duplication_node_of_objects,
                          pair_of_objects_yielding_duplication_node_distance);

      if (leftId == otherLeftId) {
        dup_node_query << "INSERT INTO "
          << "duplication_distance_of_duplication_node "
          << "(tree_id, "
          << "left_id, "
          << "duplication_distance, "
          << "first_duplicated_" << object_name << "_id, "
          << "second_duplicated_" << object_name << "_id, "
          << "left_id_where_duplicated, "
          << "distance_from_node_where_duplicated_to_first, "
          << "distance_from_node_where_duplicated_to_second) "
          << " VALUES ("
          << tree_id << ", "
          << leftId << ", "
          << duplication_distance << ", "
          << lesserObj << ", "
          << greaterObj << ", "
          << otherLeftId << ", "
          << dupInfoPtr->getDistanceToLesserObj() << ", "
          << dupInfoPtr->getDistanceToGreaterObj() << ")";
      } else {
        dup_node_query << "INSERT INTO "
          << "duplication_distance_of_duplication_node "
          << "(tree_id, "
          << "left_id, "
          << "duplication_distance, "
          << "first_duplicated_" << object_name << "_id, "
          << "second_duplicated_" << object_name << "_id, "
          << "left_id_where_duplicated) "
          << " VALUES ("
          << tree_id << ", "
          << leftId << ", "
          << duplication_distance << ", "
          << lesserObj << ", "
          << greaterObj << ", "
          << otherLeftId << ")";
      }
      T.exec(dup_node_query.str());
    }
    for (int_objobj_iter = pair_of_nearest_objects_to_maximal_node.begin();
        int_objobj_iter != pair_of_nearest_objects_to_maximal_node.end();
        ++int_objobj_iter) {
      leftId = idMap.find(int_objobj_iter->first)->second.first;
      lesserObj = int_objobj_iter->second->getLesserObj();
      greaterObj = int_objobj_iter->second->getGreaterObj();
      dupInfoPtr = duplication_node_of_objects.find(lesserObj)->second->find(
                    greaterObj)->second;
      otherLeftId 
        = idMap.find(dupInfoPtr->getDuplicationNodeId())->second.first;
      ostringstream max_node_query;
      max_node_query << "INSERT INTO "
        << "duplication_distance_of_maximal_node "
        << "(tree_id, "
        << "left_id, "
        << "first_maximizing_" << object_name << "_id, "
        << "second_maximizing_" << object_name << "_id, "
        << "left_id_where_duplicated) "
        << " VALUES ("
        << tree_id << ", "
        << leftId << ", "
        << lesserObj << ", "
        << greaterObj << ", "
        << otherLeftId << ")";
      T.exec(max_node_query.str());
    }
  }
  T.commit();
  return alreadyInDatabase;
}
