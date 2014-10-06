// get_orthologs_from_db.h
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

#ifndef GET_ORTHOLOGS_FROM_DB_H
#define GET_ORTHOLOGS_FROM_DB_H

#include <Phyl/Tree.h>
#include <map>
#include <string>
#include "modified_preorder_tree_traversal.h"
#include "find_orthologs_in_tree.h"
#include "orthologs_common.h"
#ifdef USING_POSTGRES
#include <pqxx/connection>
#else
#include <mysql++/connection.h>
#endif

//template<class ObjectT, ObjectT nullObjectValue,
//        class AttributeT, AttributeT nullAttributeValue>
typedef DBIdentifierT AttributeT;
typedef DBIdentifierT ObjectT;
void getOrthologsFromDB(const bpp::Tree &tree,
            const map<int, int> &nodeIdOfLeftIdMap,
            const int tree_id,
            const string &attribute_name,
            const string &object_name,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
            map<AttributeT , ObjectT > 
                &unique_object_with_attribute,
            map<int, map<AttributeT , TreeDistanceInfo<ObjectT, 
                                                        nullObjectValue> *> *>
                &distance_from_object_with_attribute_to_node,
            map<int, set<AttributeT > *> 
                &attributes_with_multiple_objects_for_which_node_is_maximal);
#endif
