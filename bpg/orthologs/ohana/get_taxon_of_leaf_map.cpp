// get_taxon_of_leaf_map.cpp
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

#include <sstream>
using std::ostringstream;

#include "get_taxon_of_leaf_map.h"

void getTaxonOfLeafMap(
    const map<int, int> &protein_sequence_of_leaf,
    map<int, int> &taxon_of_leaf, pqxx::transaction<>& xaction) 
{
    taxon_of_leaf.clear();
    
    int taxonId;
    map<int, int>::const_iterator iter;
    
    for(iter = protein_sequence_of_leaf.begin(); iter != protein_sequence_of_leaf.end(); ++iter) 
    {
        ostringstream query;

        query << "SELECT taxon_id FROM sequence_header WHERE id=" << (*iter).second;

        pqxx::result res = xaction.exec(query.str());
        
        if(res.size() > 0) 
        {
            string taxon_str = res[0]["taxon_id"].c_str();
            taxonId = (taxon_str == ""  ? 0: atoi(taxon_str.c_str()));
        } 
        else 
            taxonId = 0;
      
        taxon_of_leaf[(*iter).first] = taxonId;
    }

}
