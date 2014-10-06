// write_phyloxml_file.cpp
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

#include "write_phyloxml_file.h"
#include <fstream>
#include <iostream>

void indent(ofstream &outf, int indentation_level) {
  for (int i = 0; i < 2 * indentation_level; ++i) {
    outf << " ";
  }
}
  
void write_node_to_phyloxml(const bpp::Tree &tree,
    const LeftRightIdsOfNodeIdMap &idMap,
    const map<ObjectT, int> &leaf_of_object,
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *> &pair_of_nearest_objects_to_maximal_node,
    const map<int,  double> &greatest_distance_of_maximal_descendant,
    const map<int, AttributeT > &attribute_of_node,
    const vector<string> &attribute_names,
    vector<pair<int, int> > &duplicationNodeLeftIdOfMaximalNodeLeftId,
    vector<pair<int, int> > &descendantLeftIdYieldingDuplicationNodeDistance,
    vector<pair<int, pair<int, int> > > &objectLeftIdsOfMaximalNodeLeftId,
    vector<pair<int, pair<int, int> > > &objectLeftIdsOfDuplicationNodeLeftId,
    ofstream &outf, int nodeId, int &indentation_level) {
  int leftId = 0;
  indent(outf, indentation_level++);
  leftId = idMap.find(nodeId)->second.first;
  outf << "<clade id_source=\"left_" << leftId << "\">" << endl;
  if (tree.isLeaf(nodeId)) {
    indent(outf, indentation_level);
    outf << "<name>" << tree.getNodeName(nodeId) << "</name>" << endl;
  }
  if (!tree.isRoot(nodeId)) {
    indent(outf, indentation_level);
    outf << "<branch_length>" << tree.getDistanceToFather(nodeId);
    outf << "</branch_length>" << endl;
  }
  map<int, DBIdentifierT>::const_iterator attr_iter;
  attr_iter = attribute_of_node.find(nodeId);
  if (tree.isLeaf(nodeId) && attr_iter != attribute_of_node.end()) {
    indent(outf, indentation_level++);
    outf << "<taxonomy>" << endl;
    indent(outf, indentation_level);
    outf << "<code>" << attribute_names[attr_iter->second] << "</code>" << endl;
    indent(outf, --indentation_level);
    outf << "</taxonomy>" << endl;
  }
  map<int, ObjObjInfo *>::const_iterator dup_node_iter;
  dup_node_iter = pair_of_objects_yielding_duplication_node_distance.find(nodeId);
  if (dup_node_iter != pair_of_objects_yielding_duplication_node_distance.end()) {
    indent(outf, indentation_level++);
    outf << "<events>" << endl;
    indent(outf, indentation_level);
    outf << "<duplications>1</duplications>" << endl;
    indent(outf, --indentation_level);
    outf << "</events>" << endl;
    indent(outf, indentation_level);
    outf << "<property datatype=\"xsd:double\" ref=\"DUP:duplication_distance\" "
         << "applies_to=\"clade\">";
    outf << getDuplicationNodeDistance(nodeId, duplication_node_of_objects,
              pair_of_objects_yielding_duplication_node_distance);
    outf << "</property>" << endl;
    const ObjectT &maximizingObj1 = dup_node_iter->second->getLesserObj();
    const ObjectT &maximizingObj2 = dup_node_iter->second->getGreaterObj();
    DupInfo *dupInfoPtr = duplication_node_of_objects.find(
                            maximizingObj1)->second->find(
                              maximizingObj2)->second;
    int descDupNodeId = dupInfoPtr->getDuplicationNodeId();
    int descDupNodeLeftId = idMap.find(descDupNodeId)->second.first;
    descendantLeftIdYieldingDuplicationNodeDistance.push_back(make_pair(leftId,
                                                            descDupNodeLeftId));
    if (descDupNodeId == nodeId) {
      int obj1LeftId 
        = idMap.find(leaf_of_object.find(maximizingObj1)->second)->second.first;
      int obj2LeftId
        = idMap.find(leaf_of_object.find(maximizingObj2)->second)->second.first;
      objectLeftIdsOfDuplicationNodeLeftId.push_back(make_pair(leftId,
                                          make_pair(obj1LeftId, obj2LeftId)));
      indent(outf, indentation_level);
      outf  << "<property datatype=\"xsd:double\" "
            << "ref=\"DUP:distance_to_first_duplicated_leaf\" "
            << "applies_to=\"clade\">";
      outf << dupInfoPtr->getDistanceToLesserObj();
      outf << "</property>" << endl;
      indent(outf, indentation_level);
      outf  << "<property datatype=\"xsd:double\" "
            << "ref=\"DUP:distance_to_second_duplicated_leaf\" "
            << "applies_to=\"clade\">";
      outf << dupInfoPtr->getDistanceToGreaterObj();
      outf << "</property>" << endl;
    }
  }
  map<int, ObjObjInfo *>::const_iterator max_node_iter;
  max_node_iter = pair_of_nearest_objects_to_maximal_node.find(nodeId);
  if (max_node_iter != pair_of_nearest_objects_to_maximal_node.end()) {
    map<int, double>::const_iterator max_desc_dist_iter;
    max_desc_dist_iter = greatest_distance_of_maximal_descendant.find(nodeId);
    const ObjectT &maximizingObj1 = max_node_iter->second->getLesserObj();
    const ObjectT &maximizingObj2 = max_node_iter->second->getGreaterObj();
    int dupNodeId = duplication_node_of_objects.find(
                            maximizingObj1)->second->find(
                              maximizingObj2)->second->getDuplicationNodeId();
    double min_threshold, max_threshold;
    min_threshold = max(max_desc_dist_iter->second, 0.0);
    max_threshold = getDuplicationNodeDistance(dupNodeId, 
                            duplication_node_of_objects,
                            pair_of_objects_yielding_duplication_node_distance);
    if (min_threshold < max_threshold) {
      indent(outf, indentation_level);
      outf << "<property datatype=\"xsd:double\" ref=\"PHOG:min_threshold\" "
              "applies_to=\"clade\">";
      outf << max(max_desc_dist_iter->second, 0.0);
      outf << "</property>" << endl;
      indent(outf, indentation_level);
      outf << "<property datatype=\"xsd:double\" ref=\"PHOG:max_threshold\" "
              "applies_to=\"clade\">";
      outf << getDuplicationNodeDistance(dupNodeId, duplication_node_of_objects,
                pair_of_objects_yielding_duplication_node_distance);
      outf << "</property>" << endl;
      int dupNodeLeftId = idMap.find(dupNodeId)->second.first;
      duplicationNodeLeftIdOfMaximalNodeLeftId.push_back(make_pair(leftId, 
                                                                  dupNodeLeftId));
      int obj1LeftId 
        = idMap.find(leaf_of_object.find(maximizingObj1)->second)->second.first;
      int obj2LeftId
        = idMap.find(leaf_of_object.find(maximizingObj2)->second)->second.first;
      objectLeftIdsOfMaximalNodeLeftId.push_back(make_pair(leftId,
                                            make_pair(obj1LeftId, obj2LeftId)));
    }
  }
  const vector<int> &children = tree.getSonsId(nodeId);
  for (int j = 0; j < children.size(); ++j) {
    write_node_to_phyloxml(tree, idMap, leaf_of_object, 
      duplication_node_of_objects,
      pair_of_objects_yielding_duplication_node_distance,
      pair_of_nearest_objects_to_maximal_node,
      greatest_distance_of_maximal_descendant,
      attribute_of_node, attribute_names, 
      duplicationNodeLeftIdOfMaximalNodeLeftId,
      descendantLeftIdYieldingDuplicationNodeDistance,
      objectLeftIdsOfMaximalNodeLeftId,
      objectLeftIdsOfDuplicationNodeLeftId,
      outf, children[j], indentation_level);
  }
  indent(outf, --indentation_level);
  outf << "</clade>" << endl;
}

void write_PhyloXML_File(const bpp::Tree &tree,
    const LeftRightIdsOfNodeIdMap &idMap,
    const map<ObjectT, int> &leaf_of_object,
    const map<ObjectT, map<ObjectT, DupInfo *> *> &duplication_node_of_objects,
    const map<int, ObjObjInfo *>
      &pair_of_objects_yielding_duplication_node_distance,
    const map<int, ObjObjInfo *> &pair_of_nearest_objects_to_maximal_node,
    const map<int,  double> &greatest_distance_of_maximal_descendant,
    const map<int, AttributeT > &attribute_of_node,
    const vector<string> &attribute_names,
    const string &xmlfilepath) {
  ofstream outf(xmlfilepath.c_str(), ios::out);
  outf << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
  outf << "<phyloxml "
            << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
            << "xmlns=\"http://www.phyloxml.org\" "
            << "xsi:schemaLocation=\"http://www.phyloxml.org "
            << "http://www.phyloxml.org/1.10/phyloxml.xsd\">"
            << endl;
  outf << "  <phylogeny rooted=\"true\">" << endl;
  int indentation_level = 2;
  vector<pair<int, int> > duplicationNodeLeftIdOfMaximalNodeLeftId;
  vector<pair<int, pair<int, int> > > objectLeftIdsOfMaximalNodeLeftId;
  vector<pair<int, int> > descendantLeftIdYieldingDuplicationNodeDistance;
  vector<pair<int, pair<int, int> > > objectLeftIdsOfDuplicationNodeLeftId;
  write_node_to_phyloxml(tree, idMap, leaf_of_object, duplication_node_of_objects,
                        pair_of_objects_yielding_duplication_node_distance,
                        pair_of_nearest_objects_to_maximal_node,
                        greatest_distance_of_maximal_descendant,
                        attribute_of_node, attribute_names, 
                        duplicationNodeLeftIdOfMaximalNodeLeftId,
                        descendantLeftIdYieldingDuplicationNodeDistance,
                        objectLeftIdsOfMaximalNodeLeftId,
                        objectLeftIdsOfDuplicationNodeLeftId,
                        outf, tree.getRootId(), indentation_level);
  for (int m = 0; m < duplicationNodeLeftIdOfMaximalNodeLeftId.size(); ++m) {
    outf  << "  <clade_relation id_ref_0=\"left_"
          << duplicationNodeLeftIdOfMaximalNodeLeftId[m].first << "\" "
          << "id_ref_1=\"left_"
          << duplicationNodeLeftIdOfMaximalNodeLeftId[m].second << "\" "
          << "type=\"duplication_node_of_PHOG\" />" << endl;
  }
  for (int m = 0; m < objectLeftIdsOfMaximalNodeLeftId.size(); ++m) {
    outf  << "  <clade_relation id_ref_0=\"left_"
          << objectLeftIdsOfMaximalNodeLeftId[m].first << "\" "
          << "id_ref_1=\"left_"
          << objectLeftIdsOfMaximalNodeLeftId[m].second.first << "\" "
          << "type=\"first_duplicated_leaf_for_PHOG\" />" << endl;
    outf  << "  <clade_relation id_ref_0=\"left_"
          << objectLeftIdsOfMaximalNodeLeftId[m].first << "\" "
          << "id_ref_1=\"left_"
          << objectLeftIdsOfMaximalNodeLeftId[m].second.second << "\" "
          << "type=\"second_duplicated_leaf_for_PHOG\" />" << endl;
  }
  for (int d = 0; d < descendantLeftIdYieldingDuplicationNodeDistance.size();
      ++d) {
    outf << "  <clade_relation id_ref_0=\"left_"
          << descendantLeftIdYieldingDuplicationNodeDistance[d].first << "\" "
          << "id_ref_1=\"left_"
          << descendantLeftIdYieldingDuplicationNodeDistance[d].second << "\" "
          << "type=\"descendant_duplication_node_maximizing_distance\" />"
          << endl;
  }
  for (int d = 0; d < objectLeftIdsOfDuplicationNodeLeftId.size(); ++d) {
    outf  << "  <clade_relation id_ref_0=\"left_"
          << objectLeftIdsOfDuplicationNodeLeftId[d].first << "\" "
          << "id_ref_1=\"left_"
          << objectLeftIdsOfDuplicationNodeLeftId[d].second.first << "\" "
          << "type=\"first_duplicated_leaf_for_duplication_node\" />" << endl;
    outf  << "  <clade_relation id_ref_0=\"left_"
          << objectLeftIdsOfDuplicationNodeLeftId[d].first << "\" "
          << "id_ref_1=\"left_"
          << objectLeftIdsOfDuplicationNodeLeftId[d].second.second << "\" "
          << "type=\"second_duplicated_leaf_for_duplication_node\" />" << endl;
  }
  outf << "  </phylogeny>" << endl;
  outf << "</phyloxml>" << endl;
  outf.close();
}
