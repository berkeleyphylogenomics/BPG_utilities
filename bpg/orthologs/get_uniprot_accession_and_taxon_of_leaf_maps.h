// get_uniprot_accession_and_taxon_of_leaf_maps.h
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

#ifndef GET_UNIPROT_ACCESSION_AND_TAXON_OF_LEAF_MAPS_H
#define GET_UNIPROT_ACCESSION_AND_TAXON_OF_LEAF_MAPS_H

#include <map>
#include <string>
#include "orthologs_common.h"
#ifdef USING_POSTGRES
#include <pqxx/connection>
#else
#include <mysql++/connection.h>
#endif
using namespace std;

extern void getUniProtAccessionAndTaxonOfLeafMaps(
    const map<int, DBIdentifierT> &sequence_header_of_leaf,
#ifdef USING_POSTGRES
										pqxx::connection &conn,
#else
										mysqlpp::Connection &conn,
#endif
    map<int, string> &uniprot_accession_of_leaf,
    map<int, DBIdentifierT> &taxon_of_leaf);

#endif
