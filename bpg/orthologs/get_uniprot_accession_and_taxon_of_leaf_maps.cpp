// get_uniprot_accession_and_taxon_of_leaf_maps.cpp
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

#include "get_uniprot_accession_and_taxon_of_leaf_maps.h"
#ifdef USING_POSTGRES
#include <pqxx/transactor>
#include <pqxx/result>
#else
#include <mysql++/query.h>
#include <mysql++/result.h>
#endif
#include <iostream>
#include <string>
#include <sstream>
using std::ostringstream;

void getUniProtAccessionAndTaxonOfLeafMaps(
    const map<int, DBIdentifierT> &sequence_header_of_leaf,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
    map<int, string> &uniprot_accession_of_leaf,
    map<int, DBIdentifierT> &taxon_of_leaf) {
  uniprot_accession_of_leaf.clear();
  taxon_of_leaf.clear();
#ifdef USING_POSTGRES
    pqxx::work T(conn, "getTaxonOfLeafMapTransaction");
    pqxx::result R;
#endif
    ostringstream sequence_header_ids_clause;
    sequence_header_ids_clause << "(";
    map<DBIdentifierT, pair<string, DBIdentifierT> > 
        uniprot_accession_and_taxon_of_sequence_header_id;
    map<int, DBIdentifierT>::const_iterator iter;
    iter = sequence_header_of_leaf.begin();
    if (iter != sequence_header_of_leaf.end()) {
      sequence_header_ids_clause << (*iter).second;
      ++iter;
    }
    while (iter != sequence_header_of_leaf.end()) {
      sequence_header_ids_clause << "," << (*iter).second;
      ++iter;
    }
    sequence_header_ids_clause << ")";
#ifdef USING_POSTGRES
    ostringstream uniprot_query;
    uniprot_query << "SELECT sequence_header.id AS sequence_header_id, "
                  << "accession, uniprot.taxon_id "
                  << "FROM sequence_header, uniprot "
                  << "WHERE sequence_header.id IN "
                  << sequence_header_ids_clause.str()
                  << " AND uniprot_id = uniprot.id";
#else
    mysqlpp::Query uniprot_query = conn.query();
    uniprot_query << "SELECT protein_sequence.id AS sequence_header_id, "
                  << "uniprot_accession, ncbi_taxid "
                  << "FROM protein_sequence, genbank_uniprot "
                  << "WHERE protein_sequence.id IN "
                  << sequence_header_ids_clause.str()
                  << " AND genbank_uniprot_id = genbank_uniprot.id";
#endif
#ifdef USING_POSTGRES
    R = T.exec(uniprot_query.str());
    int taxon;
    DBIdentifierT sequence_header_id;
    string uniprot_accession;
    for (pqxx::result::const_iterator i = R.begin(); i != R.end(); ++i) {
      (*i)[0].to(sequence_header_id);
      (*i)[1].to(uniprot_accession);
      if ((*i)[2].is_null()) {
        taxon = 0;
      } else {
        (*i)[2].to(taxon);
      }
      uniprot_accession_and_taxon_of_sequence_header_id[sequence_header_id]
        = make_pair(uniprot_accession, taxon);
    }
#else
    mysqlpp::ResUse uniprot_res = uniprot_query.use();
    mysqlpp::Row uniprot_row;
    while (uniprot_row = uniprot_res.fetch_row()) {
      string taxon_str = uniprot_row["ncbi_taxid"].get_string();
      int taxon = (taxon_str == "NULL"  ? 0
                                        : int(uniprot_row["ncbi_taxid"]));
      uniprot_accession_and_taxon_of_sequence_header_id[
                                      int(uniprot_row["sequence_header_id"])]
        = make_pair(uniprot_row["uniprot_accession"].get_string(), taxon);
    }
#endif
    for(iter = sequence_header_of_leaf.begin(); 
        iter != sequence_header_of_leaf.end(); ++iter) {
      uniprot_accession_of_leaf[(*iter).first] 
          = uniprot_accession_and_taxon_of_sequence_header_id[
                                                          (*iter).second].first;
      taxon_of_leaf[(*iter).first] 
          = uniprot_accession_and_taxon_of_sequence_header_id[
                                                          (*iter).second].second;
    }
}
