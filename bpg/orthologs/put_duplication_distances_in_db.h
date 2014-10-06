// put_duplication_distances_in_db.h
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

#ifndef PUT_DUPLICATION_DISTANCES_IN_DB_H
#define PUT_DUPLICATION_DISTANCES_IN_DB_H

#include <Phyl/Tree.h>
#include <map>
#include <string>
#include <set>
#include "modified_preorder_tree_traversal.h"
#include "find_duplication_distances.h"
#include "orthologs_common.h"
#ifdef USING_POSTGRES
#include <pqxx/connection>
#else
#include <mysql++/connection.h>
#endif
using namespace std;

//template<class ObjectT, ObjectT nullObjectValue>
typedef DBIdentifierT ObjectT;
bool putDuplicationDistancesInDB(const bpp::Tree &tree,
    const LeftRightIdsOfNodeIdMap &idMap,
    const int tree_id,
    const string object_name,
#ifdef USING_POSTGRES
    pqxx::connection &conn,
#else
    mysqlpp::Connection &conn,
#endif
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *> &pair_of_nearest_objects_to_maximal_node,
    const map<int,  double> &greatest_distance_of_maximal_descendant);
#endif
