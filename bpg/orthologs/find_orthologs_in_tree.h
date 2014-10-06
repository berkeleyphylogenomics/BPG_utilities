// find_orthologs_in_tree.h
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

#ifndef FIND_ORTHOLOGS_IN_TREE_H
#define FIND_ORTHOLOGS_IN_TREE_H

#include <Phyl/Tree.h>
#include <map>
#include <set>
#include "orthologs_common.h"

using namespace std;

template<class ObjectT, ObjectT nullObjectValue>
class TreeDistanceInfo {
  private:
    ObjectT _nearestObjectWithAttributeValue;
    ObjectT _nextNearestObjectWithAttributeValue;
    double _distanceToNearestObjectWithAttributeValue;
    double _distanceToNextNearestObjectWIthAttributeValue;
    bool _nearestObjectIsDescendant;
    bool _isObjectPure;
  public:
    TreeDistanceInfo() : _nearestObjectWithAttributeValue(nullObjectValue),
                        _nextNearestObjectWithAttributeValue(nullObjectValue),
                        _distanceToNearestObjectWithAttributeValue(-1.0),
                        _distanceToNextNearestObjectWIthAttributeValue(-1.0),
                        _nearestObjectIsDescendant(false),
                        _isObjectPure(false)
      {};
    virtual ~TreeDistanceInfo() {};
    const ObjectT getNearestObjectWithAttributeValue() const {
      return _nearestObjectWithAttributeValue;
    };
    const ObjectT getNextNearestObjectWithAttributeValue() const {
      return _nextNearestObjectWithAttributeValue;
    };
    double getDistanceToNearestObjectWithAttributeValue() const {
      return _distanceToNearestObjectWithAttributeValue;
    };
    double getDistanceToNextNearestObjectWithAttributeValue() const {
      return _distanceToNextNearestObjectWIthAttributeValue;
    };
    bool isNearestObjectDescendant() const {
      return _nearestObjectIsDescendant;
    };
    bool isObjectPure() const {
      return _isObjectPure;
    };
    void setNearestObjectWithAttributeValue(const ObjectT object) {
      _nearestObjectWithAttributeValue = object;
    };
    void setNextNearestObjectWithAttributeValue(const ObjectT object) {
      _nextNearestObjectWithAttributeValue = object;
    };
    void setDistanceToNearestObjectWithAttributeValue(double distance) {
      _distanceToNearestObjectWithAttributeValue = distance;
    };
    void setDistanceToNextNearestObjectWithAttributeValue(double distance) {
      _distanceToNextNearestObjectWIthAttributeValue = distance;
    };
    void setNearestObjectIsDescendant(bool nearestObjectIsDescendant) {
      _nearestObjectIsDescendant = nearestObjectIsDescendant;
    };
    void setIsObjectPure(bool isObjectPure) {
      _isObjectPure = isObjectPure;
    };
};

//template<class ObjectT, ObjectT nullObjectValue, class AttributeT, 
//          AttributeT nullAttributeValue>
typedef DBIdentifierT AttributeT;
extern DBIdentifierT nullAttributeValue;
typedef DBIdentifierT ObjectT;
const DBIdentifierT nullObjectValue = 0;
void find_orthologs_in_tree(const bpp::Tree &tree, 
        const vector<int> &breadth_first_visit_order,
        const map<int, ObjectT> &object_of_leaf,
        const map<int, AttributeT> &attribute_of_node,
        map<AttributeT, ObjectT> &unique_object_with_attribute,
        set<AttributeT> &attributes_with_multiple_objects,
        map<int, map<AttributeT, TreeDistanceInfo<ObjectT, 
                                                    nullObjectValue> *> *>
          &distance_from_object_with_attribute_to_node,
          bool check_attributes = true,
          bool seek_longest = false);
#endif
