#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "modified_preorder_tree_traversal.h"
#include "get_sequence_header_of_leaf_map.h"
#include "get_taxon_of_leaf_map.h"
#include "get_breadth_first_visit_order.h"
#include "get_longest_distance_in_tree.h"
#include "find_inparalogs_in_tree.h"
#include "put_inparalogs_in_db.h"
#include "find_orthologs_in_tree.h"
#include "find_proximal_subtrees.h"
#include "find_duplication_distances.h"
#include "put_orthologs_in_db.h"
#include "put_superorthologous_nodes_in_db.h"
#include "put_duplication_distances_in_db.h"

#ifdef USING_POSTGRES
#include <unistd.h>
#include <sys/param.h>
#include <fstream>
#include <mysql++/common.h>
#include "get_tree_node_and_sequence_header_ids.h"
#include <pqxx/transactor>
#include <pqxx/result>
#else
#include <mysql++/connection.h>
#include <mysql++/query.h>
#include <mysql++/result.h>
#include "put_left_right_ids_in_db.h"
#endif

int main(int ac, char* av[])
{
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("book", po::value<string>(), "book within which to find orthologs")
      ("method", po::value<string>(), "construction method of tree to use")
      ("database", po::value<string>(), "database to connect to")
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    if (vm.count("book")) {
#ifdef USING_POSTGRES
        /* DB access parameters */
        string host = "db";
        string database = "pfacts003_test";
        if (vm.count("database")) {
            database = vm["database"].as<string>();
        } 
        string dbaccess = "host=" + host 
                        + " dbname=" + database;
        /* Read the credentials file to retrieve the password */
        ifstream inFile("/clusterfs/ohana/software/db/bpg_credentials.cfg");
        bool found_password = false;
        bool found_user = false;
        string password;
        char buffer[256];
        while (!found_password and !inFile.eof()) {
          inFile.getline(buffer, 255);
          if (!found_user) {
            if (strcmp(buffer, "[bpg_user]") == 0) {
              found_user = true;
            }
          } else {
            password = string(buffer).substr(10);
            found_password = true;
          }
        }
        dbaccess += " user=bpg_user password=" + password;
        pqxx::connection conn(dbaccess);
        boost::regex bpg_patt("(bpg[0-9][0-9][0-9][0-9])[0-9][0-9][0-9]");
        boost::smatch book_matches;
        string book = vm["book"].as<string>();
        string book_dir = "/clusterfs/ohana/bpg/pfacts/";
        if (boost::regex_match(book, book_matches, bpg_patt)) {
          book_dir += book.substr(0,4) + "/" + book_matches[1] + "/" 
                      + book + "/";
        } else {
          cerr << "Invalid book specification " << book << endl;
          return 1;
        }
#else
        const char * db = "pfacts002", *server = "phylodb";
        const char *user = "pfacts002", *pass = "fritzLang";
        mysqlpp::Connection conn(false);
        if (conn.connect(db, server, user, pass)) {
          boost::regex bpg_patt("(bpg[0-9][0-9][0-9])[0-9][0-9][0-9]");
          boost::smatch book_matches;
          string book = vm["book"].as<string>();
          string book_dir = "/home/bpg/pfacts/";

          mysqlpp::Query gathering_method_query = conn.query();
          gathering_method_query << "SELECT best_gathering_method_sw2 "
                                  << "AS gathering_method "
                                  << "FROM book "
                                  << "WHERE scopid_or_bpgid = '" << book << "'";
          mysqlpp::ResUse gathering_method_res = gathering_method_query.use();
          mysqlpp::Row gathering_method_row = gathering_method_res.fetch_row();
          if (gathering_method_row == 0) {
            cout << "Book " << book << " not found in database." << endl;
            exit(1);
          }
          string gathering_method 
            = gathering_method_row["gathering_method"].get_string();
          if (gathering_method == "NULL" || gathering_method == "") {
            gathering_method = "user";
          }
          while (gathering_method_res.fetch_row()) {
            ;
          }
              
          if (boost::regex_match(book, book_matches, bpg_patt)) {
            book_dir += book_matches[1] + "/" + book + "/" + gathering_method 
                                        + "/";
          } else {
            book_dir += book + "/" + gathering_method + "/";
          }
#endif
          string method;
          if (vm.count("method")) {
            method = vm["method"].as<string>();
          } else {
#ifdef USING_POSTGRES
            pqxx::work MT(conn, "findCanonicalTreeMethodTransaction");
            pqxx::result MR;
            ostringstream method_query;
            method_query  << "SELECT method "
                          << "FROM family, tree "
                          << "WHERE family.id = " << book.substr(3)
                          << "AND family_id = family.id "
                          << "AND canonical_tree_id = tree.id";
            MR = MT.exec(method_query.str());
            if (MR.size() == 0) {
              cerr << "Canonical tree method for " << book 
                    << " not found in database." << endl;
              return 1;
            }
            // canonical_tree_id is constrained to be not null
            MR[0][0].to(method);
#else
            method = "nj";
#endif
          }
          bpp::Newick * newickReader = new bpp::Newick(false);
          bpp::Tree * tree;
          tree = newickReader->read(book_dir + book + "." + method);
          int max_right_id;
          LeftRightIdsOfNodeIdMap idMap;
          map<int, int> nodeIdOfLeftIdMap;
          map<int, int> levelOfNodeIdMap;
          cout << "Mapping node ids to left & right ids" << endl;
          max_right_id = findLeftRightIds(*tree, idMap, nodeIdOfLeftIdMap,
                                          levelOfNodeIdMap,
                                          tree->getRootId(), 1, 0);
          map<int, DBIdentifierT> db_id_of_left_id;
          map<int, DBIdentifierT> sequence_header_of_leaf;
          map<DBIdentifierT, int> leaf_of_sequence_header;
#ifdef USING_POSTGRES
          pqxx::work T(conn, "findOrthologsFindTreeTransaction");
          pqxx::result R;
          DBIdentifierT tree_id;
          {
            ostringstream query;
            query << "SELECT id "
                  << "FROM tree "
                  << "WHERE family_id = " << book.substr(3)
                  << " AND method = '" << method << "'";
            R = T.exec(query.str());
            if (R.size() == 0) {
              cerr  << "No " << method << " tree for family " << book
                    << " found in database." << endl;
              return 1;
            }
            R[0][0].to(tree_id);
          }
          ostringstream query;
          query << "SELECT reciprocal_longest_distance "
                << "FROM tree "
                << "WHERE id = " << tree_id;
          R = T.exec(query.str());
          T.commit();
          double reciprocal_longest_distance;
          if (R[0][0].is_null()) {
            reciprocal_longest_distance = -1.0;
          } else {
            R[0][0].to(reciprocal_longest_distance);
          }
          if (!getTreeNodeAndSequenceHeaderIds(*tree, idMap, 
                              nodeIdOfLeftIdMap, 
                              levelOfNodeIdMap,
                              tree_id, 
                              conn,
                              db_id_of_left_id,
                              sequence_header_of_leaf,
                              leaf_of_sequence_header)) {
            cerr  << "An error occurred while retrieving tree nodes from db."
                  << endl;
            return 1;
          }
#else
          string alignment_id;
          mysqlpp::Query alignment_query = conn.query();
          alignment_query << "SELECT alignment.id AS alignment_id "
                          << "FROM alignment, book "
                          << "WHERE book_id = book.id "
                          << "AND subfamily_alignment_id = 0 "
                          << "AND gathering_method = '" 
                          << gathering_method << "' "
                          << "AND scopid_or_bpgid = '" << book << "'";
          mysqlpp::ResUse alignment_res = alignment_query.use();
          mysqlpp::Row alignment_row = alignment_res.fetch_row();
          alignment_id = alignment_row["alignment_id"].get_string();
          while (alignment_res.fetch_row()) {
            ;
          }
          cout << "Making sure tree is in the database" << endl;
          DBIdentifierT tree_id;
          mysqlpp::Query insert_tree_query = conn.query();
          insert_tree_query << "INSERT INTO tree "
                            << "(alignment_id, method, is_rsd_rooted) "
                            << "VALUES "
                            << "(" << alignment_id << ", 'nj', FALSE)";
          mysqlpp::Query tree_query = conn.query();
          ostringstream tree_query_str;
          tree_query_str << "SELECT id, reciprocal_longest_distance FROM tree "
                      << "WHERE alignment_id = " << alignment_id
                      << " AND method = 'nj' "
                      << "AND is_rsd_rooted = FALSE";
          tree_query << tree_query_str.str();
          mysqlpp::ResUse tree_res = tree_query.use();
          mysqlpp::Row tree_row = tree_res.fetch_row();
          if (!tree_row) {
            mysqlpp::ResNSel status = insert_tree_query.execute();
            if (conn.errnum()) {
              cerr << "Error: " << conn.error() << endl;
            }
            tree_query << tree_query_str.str();
            tree_res = tree_query.use();
            tree_row = tree_res.fetch_row();
          }
          tree_id = int(tree_row["id"]);
          double reciprocal_longest_distance 
            = double(tree_row["reciprocal_longest_distance"]);
          while (tree_row = tree_res.fetch_row()) {
            ;
          }
          cout << "Storing left & right ids & protein sequence ids in db";
          cout << endl;
          if (!putLeftRightIdsInDB(*tree, idMap, nodeIdOfLeftIdMap, 
                              levelOfNodeIdMap,
                              tree_id, alignment_id,
                              conn,
                              db_id_of_left_id,
                              sequence_header_of_leaf)) {
            cout << "Left & right ids & protein sequence ids "
                  << "were already in the db, did not store " << endl;

          }
#endif
          map<int, DBIdentifierT> species_of_node;
          cout << "Retrieving species at leaves" << endl;
          getTaxonOfLeafMap(sequence_header_of_leaf, conn, species_of_node);
          cout << "Getting breadth-first visit order" << endl;
          vector<int> nodes_to_visit;
          getBreadthFirstVisitOrder(*tree, nodes_to_visit);
          if (reciprocal_longest_distance == -1.0) {
            cout << "Reciprocal longest distance in tree was not in database, "
                  << "computing... " << endl;
            double longest_distance = get_longest_distance_in_tree(*tree,
                                                              nodes_to_visit);
            reciprocal_longest_distance = 1.0 / longest_distance;
            cout << "Storing reciprocal longest distance "
                << reciprocal_longest_distance << " in database" << endl;
#ifdef USING_POSTGRES
            ostringstream span_query;
#else
            mysqlpp::Query span_query = conn.query();
#endif
            span_query << "UPDATE tree "
                       << "SET reciprocal_longest_distance = "
                       << reciprocal_longest_distance
                       << " WHERE id = " << tree_id;
#ifdef USING_POSTGRES
            pqxx::work DT(conn, "UpdateReciprocalLongestDistanceTransaction");
            DT.exec(span_query.str());
            DT.commit();
#else
            mysqlpp::ResNSel status = span_query.execute();
            if (conn.errnum()) {
              cerr << "Error: " << conn.error() << endl;
            }
#endif
          }
          cout << "Finding inparalogs in tree" << endl;
          find_inparalogs_in_tree(*tree, nodes_to_visit, species_of_node);
          cout << "Storing inparalogs in database" << endl;
          if (!putInparalogsInDB(*tree, idMap, db_id_of_left_id,
                                "species", conn, species_of_node)) {
            cout << "Inparalogs were already in database, didn't store" << endl;
          }
          cout << "Finding orthologs in tree" << endl;
          map<DBIdentifierT, DBIdentifierT> 
            unique_sequence_header_with_species;
          set<DBIdentifierT> species_with_multiple_sequence_headers;
          map<int, map<DBIdentifierT, TreeDistanceInfo<ObjectT,
                                            nullObjectValue> *> *>
            distance_from_sequence_header_with_taxon_to_node;
          find_orthologs_in_tree(*tree, nodes_to_visit, 
                          sequence_header_of_leaf, 
                          species_of_node,
                          unique_sequence_header_with_species,
                          species_with_multiple_sequence_headers,
                          distance_from_sequence_header_with_taxon_to_node);
          cout << "Finding proximal subtrees" << endl;
          map<int, set<DBIdentifierT> *> 
            species_with_multiple_sequence_headers_for_which_node_is_maximal;
          set<int> superorthologous_nodes;
          map<DBIdentifierT, set<int> *> maximal_nodes_of_sequence_header;
          map<DBIdentifierT, map<DBIdentifierT, set<int> *> *>
              alternative_nearest_sequence_headers;
          find_proximal_subtrees(*tree, nodes_to_visit,
            species_of_node,
            unique_sequence_header_with_species,
            species_with_multiple_sequence_headers,
            distance_from_sequence_header_with_taxon_to_node,
            idMap,
            leaf_of_sequence_header,
            sequence_header_of_leaf,
            species_with_multiple_sequence_headers_for_which_node_is_maximal,
            superorthologous_nodes,
            maximal_nodes_of_sequence_header,
            alternative_nearest_sequence_headers);
          cout << "There are " << unique_sequence_header_with_species.size()
              << " species with unique protein sequences in the tree." << endl;
          cout << "There are " 
    << species_with_multiple_sequence_headers_for_which_node_is_maximal.size()
            << " species with multiple protein sequences in the tree." << endl;
          /* RSD 2011/01/31
          // This records a lot of information which is no longer relevant,
          // takes a lot of space in the database and takes a long time to
          // insert, even though it is all correct.  All the information needed
          // to actually support our ortholog queries (including justifying
          // evidence) is recorded in putDuplicationDistancesInDB().  So this is
          // commented out for now.
          cout << "Putting orthologs in database" << endl;
          if (!putOrthologsInDB(*tree, idMap, db_id_of_left_id,
              "species", "sequence_header", 
              conn,
              unique_sequence_header_with_species,
              distance_from_sequence_header_with_taxon_to_node,
          species_with_multiple_sequence_headers_for_which_node_is_maximal)) {
            cout << "Orthologs were already in database, didn't store" << endl;
          }
          */
          cout << "Putting " << superorthologous_nodes.size()
              << " superorthologous nodes in database" << endl;
          putSuperorthologousNodesInDB(*tree, idMap, db_id_of_left_id, tree_id,
              "species", "sequence_header", conn, superorthologous_nodes);
          map<DBIdentifierT, map<DBIdentifierT, DupInfo *> *>
            duplication_node_of_sequence_headers;
          map<int, ObjObjInfo *>
            pair_of_sequence_headers_yielding_duplication_distance;
          map<int, ObjObjInfo *>
            pair_of_nearest_sequence_headers_to_maximal_node;
          map<int, double>
            greatest_distance_of_maximal_descendant;
          cout << "Finding duplication distances" << endl;
          find_duplication_distances(*tree, nodes_to_visit,
            species_of_node, distance_from_sequence_header_with_taxon_to_node,
            idMap, nodeIdOfLeftIdMap, leaf_of_sequence_header, 
            sequence_header_of_leaf,
            alternative_nearest_sequence_headers,
            duplication_node_of_sequence_headers,
            pair_of_sequence_headers_yielding_duplication_distance,
            pair_of_nearest_sequence_headers_to_maximal_node,
            greatest_distance_of_maximal_descendant);
          cout << "Putting duplication distances and supporting evidence "
            << " for them in database" << endl;
          if (!putDuplicationDistancesInDB(*tree, idMap, tree_id,
                "sequence_header", conn, duplication_node_of_sequence_headers,
                pair_of_sequence_headers_yielding_duplication_distance,
                pair_of_nearest_sequence_headers_to_maximal_node,
                greatest_distance_of_maximal_descendant)) {
            cout << "Supporting evidence was already in database, "
                << "didn't store" << endl;
          }
          cout << "Freeing memory" << endl;
          map<DBIdentifierT, map<DBIdentifierT, DupInfo *> *>::const_iterator
            obj_obj_dup_iter;
          map<DBIdentifierT, DupInfo *>::const_iterator obj_dup_iter;
          for (obj_obj_dup_iter = duplication_node_of_sequence_headers.begin();
              obj_obj_dup_iter != duplication_node_of_sequence_headers.end();
              ++obj_obj_dup_iter) {
            for (obj_dup_iter = obj_obj_dup_iter->second->begin();
                obj_dup_iter != obj_obj_dup_iter->second->end();
                ++obj_dup_iter) {
              delete obj_dup_iter->second;
            }
            delete obj_obj_dup_iter->second;
          }
          map<int, ObjObjInfo *>::const_iterator int_objobj_iter;
          for (int_objobj_iter
              = pair_of_sequence_headers_yielding_duplication_distance.begin();
              int_objobj_iter
              != pair_of_sequence_headers_yielding_duplication_distance.end();
              ++int_objobj_iter) {
            delete int_objobj_iter->second;
          }
          for (int_objobj_iter 
                = pair_of_nearest_sequence_headers_to_maximal_node.begin();
              int_objobj_iter
                != pair_of_nearest_sequence_headers_to_maximal_node.end();
              ++int_objobj_iter) {
            delete int_objobj_iter->second;
          }
          map<DBIdentifierT, set<int> *>::const_iterator obj_nodes_iter;
          for (obj_nodes_iter = maximal_nodes_of_sequence_header.begin();
              obj_nodes_iter != maximal_nodes_of_sequence_header.end();
              ++obj_nodes_iter) {
            delete (*obj_nodes_iter).second;
          }
          map<DBIdentifierT, map<DBIdentifierT, set<int> *> *>::const_iterator
            obj_objs_nodes_iter;
          for (obj_objs_nodes_iter 
                = alternative_nearest_sequence_headers.begin();
              obj_objs_nodes_iter
                != alternative_nearest_sequence_headers.end();
              ++obj_objs_nodes_iter) {
            for (obj_nodes_iter = (*obj_objs_nodes_iter).second->begin();
                obj_nodes_iter != (*obj_objs_nodes_iter).second->end();
                ++obj_nodes_iter) {
              delete (*obj_nodes_iter).second;
            }
            delete (*obj_objs_nodes_iter).second;
          }
          map<int, map<DBIdentifierT, TreeDistanceInfo<ObjectT,
                                      nullObjectValue> *> *>::const_iterator
            dist_attr_node_iter;
          map<DBIdentifierT, 
              TreeDistanceInfo<ObjectT, nullObjectValue> *>::const_iterator 
            attr_node_iter;
          for (dist_attr_node_iter 
                  = distance_from_sequence_header_with_taxon_to_node.begin();
                dist_attr_node_iter 
                  != distance_from_sequence_header_with_taxon_to_node.end();
                ++dist_attr_node_iter) {
            for (attr_node_iter = ((*dist_attr_node_iter).second)->begin();
                  attr_node_iter != ((*dist_attr_node_iter).second)->end();
                  ++attr_node_iter) {
              delete (*attr_node_iter).second;
            }
            delete (*dist_attr_node_iter).second;
          }
          cout << "Done." << endl;
          delete tree;
          delete newickReader;
#ifndef USING_POSTGRES
        }
#endif

    } else {
      cout << "Book was not set.\n";
    }
  }
  catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      return 1;
  }
  catch(...) {
      cerr << "Exception of unknown type!\n";
  }

  return 0;
}
