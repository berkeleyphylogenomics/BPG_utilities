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

#include "get_taxon_of_leaf_map.h"
#ifdef USING_POSTGRES
#include <iostream>
#include <pqxx/transactor>
#include <pqxx/result>
#else
#include <mysql++/connection.h>
#include <mysql++/query.h>
#include <mysql++/result.h>
#endif

void getTaxonOfLeafMap(
    const map<int, DBIdentifierT> &sequence_header_of_leaf,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
    map<int, DBIdentifierT> &taxon_of_leaf) {
  taxon_of_leaf.clear();
#ifdef USING_POSTGRES
    pqxx::work T(conn, "getTaxonOfLeafMapTransaction");
    pqxx::result R;
    DBIdentifierT taxonId;
    map<int, DBIdentifierT>::const_iterator iter;
    for(iter = sequence_header_of_leaf.begin(); 
        iter != sequence_header_of_leaf.end(); ++iter) {
      ostringstream query;
      query << "SELECT taxon_id "
            << "FROM sequence_header "
            << "WHERE id = "
            << (*iter).second;
      R = T.exec(query.str());
      taxonId = 0;
      if (R.size() > 0) {
        if( R[0][0].is_null()) {
          ostringstream query2;
          query2 << "SELECT uniprot.taxon_id "
                 << "FROM sequence_header, uniprot "
                 << "WHERE uniprot_id = uniprot.id "
                 << "AND sequence_header.id = "
                 << (*iter).second;
          R = T.exec(query2.str());
          if (R.size() > 0) {
            R[0][0].to(taxonId);
          }
        } else {
          R[0][0].to(taxonId);
        }
      }
      taxon_of_leaf[(*iter).first] = taxonId;
    }
#else
    mysqlpp::Query taxon_query = conn.query();
    taxon_query << "SELECT ncbi_taxid "
                << "FROM protein_sequence, genbank_uniprot "
                << "WHERE protein_sequence.id = %0q "
                << "AND genbank_uniprot_id = genbank_uniprot.id";
    taxon_query.parse();
    DBIdentifierT taxonId;
    map<int, DBIdentifierT>::const_iterator iter;
    for(iter = sequence_header_of_leaf.begin(); 
        iter != sequence_header_of_leaf.end(); ++iter) {
      mysqlpp::ResUse taxa_res = taxon_query.use((*iter).second);
      mysqlpp::Row taxon_row = taxa_res.fetch_row();
      if (taxon_row) {
        string taxon_str = taxon_row["ncbi_taxid"].get_string();
        taxonId = (taxon_str == "NULL"  ? 0
                                          : int(taxon_row["ncbi_taxid"]));
      } else {
        taxonId = 0;
      }
      taxon_of_leaf[(*iter).first] = taxonId;
      while (taxa_res.fetch_row()) {
        ;
      }
    }
#endif
}
