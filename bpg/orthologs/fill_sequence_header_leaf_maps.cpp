// fill_sequence_header_leaf_maps.cpp
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

#include "fill_sequence_header_leaf_maps.h"

using namespace std;

void fillSequenceHeaderLeafMaps(const bpp::Tree &tree,
              vector<string> &leaf_names,
              map<string, DBIdentifierT> &sequence_header_of_leaf_name,
              map<int, DBIdentifierT> &sequence_header_of_leaf,
              map<DBIdentifierT, int> &leaf_of_sequence_header) {
  leaf_names.clear();
  sequence_header_of_leaf_name.clear();
  sequence_header_of_leaf.clear();
  leaf_of_sequence_header.clear();
  const vector<int> &leaves = tree.getLeavesId();
  vector<int>::const_iterator leaf_iter;
  DBIdentifierT sequence_header;
  for (leaf_iter = leaves.begin(); leaf_iter != leaves.end(); ++leaf_iter) {
    string leaf_name = tree.getNodeName(*leaf_iter);
    sequence_header = leaf_names.size();
    leaf_names.push_back(leaf_name);
    sequence_header_of_leaf_name[leaf_name] = sequence_header;
    sequence_header_of_leaf[*leaf_iter] = sequence_header;
    leaf_of_sequence_header[sequence_header] = *leaf_iter;
  }
}
