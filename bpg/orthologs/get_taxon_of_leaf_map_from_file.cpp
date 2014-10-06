// get_taxon_of_leaf_map_from_file.cpp
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

#include "get_taxon_of_leaf_map_from_file.h"
#include <fstream>

void getTaxonOfLeafMapFromFile(const string &taxa_file, 
        const map<DBIdentifierT, int> &leaf_of_sequence_header,
        const map<string, DBIdentifierT> &sequence_header_of_leaf_name,
        map<int, DBIdentifierT> &species_of_node,
        vector<string> &taxon_names,
        map<string, DBIdentifierT> &taxon_of_taxon_name) {
  char buffer[512];
  streamsize lineLength;
  ifstream inFile(taxa_file.c_str());
  char *taxon_ptr;
  string gene_name;
  string taxon_name;
  DBIdentifierT taxon_identifier;
  map<string, DBIdentifierT>::const_iterator taxon_name_iter;
  map<string, DBIdentifierT>::const_iterator seqhdr_name_iter;
  map<DBIdentifierT, int>::const_iterator leaf_seqhdr_iter;
  while (!inFile.eof()) {
    inFile.getline(buffer, 512);
    lineLength = inFile.gcount();
    if (lineLength == 0) {
      break;
    }
    // The line should be of the form
    // <gene name>\t<taxon name>
    // Make the null-terminated string in buffer be the gene name, and the
    // null-terminated string pointed to by taxon_ptr be the taxon name.
    if (buffer[lineLength-1] == '\n') {
      buffer[lineLength-1] = '\0';
    }
    for (taxon_ptr = buffer; *taxon_ptr && *taxon_ptr != '\t'; ++taxon_ptr) 
      ;
    *taxon_ptr = '\0';
    ++taxon_ptr;
    taxon_name = taxon_ptr;
    taxon_name_iter = taxon_of_taxon_name.find(taxon_name);
    if (taxon_name_iter == taxon_of_taxon_name.end()) {
      taxon_identifier = taxon_names.size();
      taxon_names.push_back(taxon_name);
      taxon_of_taxon_name[taxon_name] = taxon_identifier;
    } else {
      taxon_identifier = taxon_name_iter->second;
    }
    gene_name = buffer;
    seqhdr_name_iter = sequence_header_of_leaf_name.find(gene_name);
    leaf_seqhdr_iter = leaf_of_sequence_header.find(seqhdr_name_iter->second);
    species_of_node[leaf_seqhdr_iter->second] = taxon_identifier;
  }
}
