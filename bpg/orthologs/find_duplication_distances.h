// find_duplication_distances.h
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

#ifndef FIND_DUPLICATION_DISTANCES_H
#define FIND_DUPLICATION_DISTANCES_H

#include <Phyl/Tree.h>
#include "orthologs_common.h"
#include "modified_preorder_tree_traversal.h"
#include "find_orthologs_in_tree.h"

//template<class AttributeT, ObjectT, ObjectT nullObjectValue>
class DupInfo {
  private:
    int _duplicationNodeId;
    double _distanceToLesserObj;
    double _distanceToGreaterObj;
  public:
    DupInfo(int duplicationNodeId, double distanceToLesserObj,
            double distanceToGreaterObj) : 
            _duplicationNodeId(duplicationNodeId),
            _distanceToLesserObj(distanceToLesserObj),
            _distanceToGreaterObj(distanceToGreaterObj)
      {};
    virtual ~DupInfo() {};
    const int getDuplicationNodeId() const {
      return _duplicationNodeId;
    };
    const double getDistanceToLesserObj() const {
      return _distanceToLesserObj;
    };
    const double getDistanceToGreaterObj() const {
      return _distanceToGreaterObj;
    };
    const double getHalfSpan() const {
      return (_distanceToLesserObj + _distanceToGreaterObj) / 2;
    };
};

class ObjObjInfo {
  private:
    ObjectT _lesserObj;
    ObjectT _greaterObj;
  public:
    ObjObjInfo() :  _lesserObj(nullObjectValue), 
                    _greaterObj(nullObjectValue)
      {};
    virtual ~ObjObjInfo() {};
    const ObjectT &getLesserObj() const {
      return _lesserObj;
    };
    const ObjectT &getGreaterObj() const {
      return _greaterObj;
    };
    void setValues(const ObjectT &lesserObj, const ObjectT &greaterObj) {
      _lesserObj = lesserObj;
      _greaterObj = greaterObj;
    };
};

double getDuplicationNodeDistance(int nodeId,
  const map<ObjectT, map<ObjectT, DupInfo *> *> 
      &duplication_node_of_objects,
  const map<int, ObjObjInfo *>
          &pair_of_objects_yielding_duplication_node_distance);

double getMaximalNodeDistance(int nodeId,
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *>
      &pair_of_nearest_objects_to_maximal_node);

void find_duplication_distances(const bpp::Tree &tree,
        const vector<int> &breadth_first_visit_order,
        const map<int, AttributeT> attribute_of_node,
        const map<int, map<AttributeT, TreeDistanceInfo<ObjectT, 
                                                    nullObjectValue> *> *>
          distance_from_object_with_attribute_to_node,
        const LeftRightIdsOfNodeIdMap &idmap,
        const map<int, int> &nodeIdOfLeftIdMap,
        const map<ObjectT, int> &leaf_of_object,
        const map<int, ObjectT> &object_of_leaf,
        const map<ObjectT, map<ObjectT, set<int> *> *> 
          &alternative_nearest_objects,
        map<ObjectT, map<ObjectT, DupInfo *> *> 
            &duplication_node_of_objects,
        map<int, ObjObjInfo *>
          &pair_of_objects_yielding_duplication_node_distance,
        map<int, ObjObjInfo *>
          &pair_of_nearest_objects_to_maximal_node,
        map<int,  double> 
          &greatest_distance_of_maximal_descendant);
#endif
