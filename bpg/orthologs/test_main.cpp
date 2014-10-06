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
#include "get_uniprot_accession_and_taxon_of_leaf_maps.h"
#include "get_breadth_first_visit_order.h"
#include "get_longest_distance_in_tree.h"
#include "find_inparalogs_in_tree.h"
#include "put_inparalogs_in_db.h"
#include "get_inparalogs_from_db.h"
#include "find_orthologs_in_tree.h"
#include "find_proximal_subtrees.h"
#include "find_duplication_distances.h"
#include "put_orthologs_in_db.h"
#include "get_orthologs_from_db.h"
#include "put_superorthologous_nodes_in_db.h"
#include "get_superorthologous_nodes_from_db.h"
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

void printNodeId(bpp::Tree &tree, int nodeId, int level, 
            LeftRightIdsOfNodeIdMap &idMap, map<int, int> &levelOfNodeIdMap,
        map<int, DBIdentifierT> &stored_sequence_header_of_leaf,
        map<int, DBIdentifierT> &retrieved_sequence_header_of_leaf,
        map<int, string> &uniprot_accession_of_leaf,
        map<DBIdentifierT, string> uniprot_accession_of_sequence_header,
        map<DBIdentifierT, string> &names_of_taxa,
        map<int, DBIdentifierT> &stored_species_of_node,
        map<int, DBIdentifierT> &retrieved_species_of_node,
        map<int, map<DBIdentifierT, TreeDistanceInfo<ObjectT,
                                          nullObjectValue> *> *>
        &stored_distance_from_sequence_header_with_taxon_to_node,
/*
        map<int, map<DBIdentifierT, TreeDistanceInfo<ObjectT,
                                          nullObjectValue> *> *>
        &retrieved_distance_from_sequence_header_with_taxon_to_node,
*/
        map<int, set<DBIdentifierT> *>
    &stored_species_with_multiple_sequence_headers_for_which_node_is_maximal,
/*
        map<int, set<DBIdentifierT> *>
&retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal,
*/
        set<int> &stored_superorthologous_nodes,
        set<int> &retrieved_superorthologous_nodes,
        map<DBIdentifierT, map<DBIdentifierT, DupInfo *> *>
          &duplication_node_of_sequence_headers,
        map<int, ObjObjInfo *>
          &pair_of_sequence_headers_yielding_duplication_distance,
        map<int, ObjObjInfo *>
          &pair_of_nearest_sequence_headers_to_maximal_node,
        map<int, double>
          &greatest_distance_of_maximal_descendant) {
  for( int i = 0; i < level; ++i) {
    cout << " ";
  }
  cout << "nodeId: " << nodeId;
  if (stored_superorthologous_nodes.find(nodeId) 
      != stored_superorthologous_nodes.end()) {
    cout << "*";
    if (retrieved_superorthologous_nodes.find(nodeId) 
        == retrieved_superorthologous_nodes.end()) {
      cout << "?";
    }
  }
  cout << " Left_Id: " << idMap[nodeId].first;
  if (tree.hasNodeName(nodeId)) {
    cout << " Name: " << tree.getNodeName(nodeId);
    if (stored_sequence_header_of_leaf[nodeId] 
        != retrieved_sequence_header_of_leaf[nodeId]) {
      cout << retrieved_sequence_header_of_leaf[nodeId];
      cout << "!=";
      cout << stored_sequence_header_of_leaf[nodeId];
    }
    cout << " Uniprot Accession: " << uniprot_accession_of_leaf[nodeId];
  }
  map<int, ObjObjInfo *>::const_iterator int_objobj_iter;
  DupInfo *dupInfoPtr;
  int_objobj_iter 
    = pair_of_sequence_headers_yielding_duplication_distance.find(nodeId);
  if (int_objobj_iter 
    != pair_of_sequence_headers_yielding_duplication_distance.end()) {
    dupInfoPtr = (*duplication_node_of_sequence_headers[
                  int_objobj_iter->second->getLesserObj()])
                  [int_objobj_iter->second->getGreaterObj()];
    cout << " Dup:" << dupInfoPtr->getHalfSpan() << " ";
    cout << int_objobj_iter->second->getLesserObj() << ":";
    cout << dupInfoPtr->getDistanceToLesserObj() << " ";
    cout << int_objobj_iter->second->getGreaterObj() << ":";
    cout << dupInfoPtr->getDistanceToGreaterObj() << " ";
  }
  if (stored_species_of_node.find(nodeId) != stored_species_of_node.end()) {
    cout << " Species: " << names_of_taxa[stored_species_of_node[nodeId]];
    if (retrieved_species_of_node[nodeId] != stored_species_of_node[nodeId]) {
      cout << stored_species_of_node[nodeId];
      cout << "!=" << retrieved_species_of_node[nodeId];
    }
  }
  set<DBIdentifierT>::const_iterator attr_iter;
  if (stored_species_with_multiple_sequence_headers_for_which_node_is_maximal.find(
        nodeId) != 
    stored_species_with_multiple_sequence_headers_for_which_node_is_maximal.end()) {
    /*
    if (retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal.find(
          nodeId) ==
      retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal.end()) {
      cout << "??";
    }
    */
    for (attr_iter = 
    stored_species_with_multiple_sequence_headers_for_which_node_is_maximal[
            nodeId]->begin();
        attr_iter !=
    stored_species_with_multiple_sequence_headers_for_which_node_is_maximal[
            nodeId]->end();
        ++attr_iter) {
      cout << " (Ortholog";
      if ((*(stored_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
            [*attr_iter]->isNearestObjectDescendant()) {
        cout << "+";
        /*
        if (!(*(retrieved_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
              [*attr_iter]->isNearestObjectDescendant()) {
          cout << "!";
        }
        */
      }
      cout << ": "
          << uniprot_accession_of_sequence_header[
            (*(stored_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
                        [*attr_iter]->getNearestObjectWithAttributeValue()];
      /*
      if ((*(stored_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
                      [*attr_iter]->getNearestObjectWithAttributeValue() !=
        (*(retrieved_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
                      [*attr_iter]->getNearestObjectWithAttributeValue()) {
        cout << "!="
        << (*(retrieved_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
                      [*attr_iter]->getNearestObjectWithAttributeValue();
      }
      */
      cout << ", Species: " << names_of_taxa[*attr_iter];
      cout << " Distance: ";
      cout << (*(stored_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
                  [*attr_iter]->getDistanceToNearestObjectWithAttributeValue();
      cout << " Next-Nearest Distance: ";
      cout << (*(stored_distance_from_sequence_header_with_taxon_to_node[
                                                                    nodeId]))
                  [*attr_iter]->getDistanceToNextNearestObjectWithAttributeValue();
      cout << ")";
    }

    int_objobj_iter 
      = pair_of_nearest_sequence_headers_to_maximal_node.find(nodeId);
    if (int_objobj_iter 
        != pair_of_nearest_sequence_headers_to_maximal_node.end()) {
      dupInfoPtr = (*duplication_node_of_sequence_headers[
                    int_objobj_iter->second->getLesserObj()])
                    [int_objobj_iter->second->getGreaterObj()];
      cout << " Duplication of Sequence Header ";
      cout << int_objobj_iter->second->getLesserObj() << ":";
      cout << dupInfoPtr->getDistanceToLesserObj();
      cout << ", Sequence Header ";
      cout << int_objobj_iter->second->getGreaterObj() << ":";
      cout << dupInfoPtr->getDistanceToGreaterObj();
      cout << " at LeftId ";
      cout << idMap.find(dupInfoPtr->getDuplicationNodeId())->second.first;
    }
  }
  map<int, double>::const_iterator node_dist_iter;
  node_dist_iter = greatest_distance_of_maximal_descendant.find(nodeId);
  if (node_dist_iter != greatest_distance_of_maximal_descendant.end()) {
    cout << " Min Threshold: " << node_dist_iter->second;
  }
  cout << endl;
  const vector<int> &children = tree.getSonsId(nodeId);
  vector<int>::const_iterator node_iter;
  for (node_iter = children.begin(); node_iter != children.end(); 
      ++node_iter) {
    printNodeId(tree, *node_iter, level + 1, idMap, levelOfNodeIdMap,
    stored_sequence_header_of_leaf, 
    retrieved_sequence_header_of_leaf, 
    uniprot_accession_of_leaf, 
    uniprot_accession_of_sequence_header,
    names_of_taxa,
    stored_species_of_node, retrieved_species_of_node,
    stored_distance_from_sequence_header_with_taxon_to_node,
/*
    retrieved_distance_from_sequence_header_with_taxon_to_node,
*/
    stored_species_with_multiple_sequence_headers_for_which_node_is_maximal, 
/*
  retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal,
*/
    stored_superorthologous_nodes,
    retrieved_superorthologous_nodes,
    duplication_node_of_sequence_headers,
    pair_of_sequence_headers_yielding_duplication_distance,
    pair_of_nearest_sequence_headers_to_maximal_node,
    greatest_distance_of_maximal_descendant);
  }
  for( int i = 0; i < level; ++i) {
    cout << " ";
  }
  cout << "nodeId: " << nodeId;
  cout << " RightId: " << idMap[nodeId].second;
  if (tree.hasNodeName(nodeId)) {
    cout << " Name: " << tree.getNodeName(nodeId);
  }
  if (level != levelOfNodeIdMap[nodeId]) {
    cout << "Stored level " << levelOfNodeIdMap[nodeId] << " != "
          << "computed level " << level;
  }
  cout << endl;
}

int main(int ac, char* av[])
{
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("book", po::value<string>(), "book within which to traverse tree")
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
          map<int, DBIdentifierT> stored_sequence_header_of_leaf,
                                        retrieved_sequence_header_of_leaf;
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
          double retrieved_reciprocal_longest_distance;
          if (R[0][0].is_null()) {
            retrieved_reciprocal_longest_distance = -1.0;
          } else {
            R[0][0].to(retrieved_reciprocal_longest_distance);
          }
          if (!getTreeNodeAndSequenceHeaderIds(*tree, idMap, 
                              nodeIdOfLeftIdMap, 
                              levelOfNodeIdMap,
                              tree_id, 
                              conn,
                              db_id_of_left_id,
                              stored_sequence_header_of_leaf,
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
          double retrieved_reciprocal_longest_distance 
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
                              stored_sequence_header_of_leaf)) {
            cout << "Left & right ids & protein sequence ids "
                  << "were already in the db, did not store " << endl;

          }
#endif
          cout << "Retrieving sequence_header_of_leaf from db" << endl;
          getSequenceHeaderOfLeafMap(tree_id, nodeIdOfLeftIdMap, 
                                    conn,
                                    retrieved_sequence_header_of_leaf);
          map<int, string> uniprot_accession_of_leaf;
          map<int, DBIdentifierT> stored_species_of_node, 
                                        retrieved_species_of_node;
          cout << "Retrieving uniprot_accession_of_leaf "
                << "and stored_species_of_leaf from db" << endl;
          getUniProtAccessionAndTaxonOfLeafMaps(
              retrieved_sequence_header_of_leaf, conn, 
              uniprot_accession_of_leaf,
              stored_species_of_node);
          map<DBIdentifierT, string> uniprot_accession_of_sequence_header;
          map<int, DBIdentifierT>::const_iterator prot_seq_leaf_iter;
          for (prot_seq_leaf_iter = stored_sequence_header_of_leaf.begin();
              prot_seq_leaf_iter != stored_sequence_header_of_leaf.end();
              ++prot_seq_leaf_iter) {
            uniprot_accession_of_sequence_header[(*prot_seq_leaf_iter).second]
              = uniprot_accession_of_leaf[(*prot_seq_leaf_iter).first];
          }
          cout << "Retrieving species names from db" << endl;
          map<DBIdentifierT, string> names_of_taxa;
          map<int, DBIdentifierT>::const_iterator species_leaf_iter;
          for (species_leaf_iter = stored_species_of_node.begin();
                species_leaf_iter != stored_species_of_node.end();
                ++species_leaf_iter) {
            ostringstream taxon_number_str;
            taxon_number_str << (*species_leaf_iter).second;
            names_of_taxa[(*species_leaf_iter).second] = taxon_number_str.str();
          }
          string taxa_clause = "(";
          map<DBIdentifierT, string>::const_iterator species_iter;
          for (species_iter = names_of_taxa.begin(); 
                species_iter != names_of_taxa.end(); ++species_iter) {
            taxa_clause += (*species_iter).second;
            taxa_clause += ",";
          }
          taxa_clause += "0) ";
          cout << "Taxa: " << taxa_clause << endl;
#ifdef USING_POSTGRES
          pqxx::work ST(conn, "GetSpeciesNamesTransaction");
          ostringstream species_name_query;
#else
          mysqlpp::Query species_name_query = conn.query();
#endif
          species_name_query << "SELECT id, scientific_name, common_name "
#ifdef USING_POSTGRES
                            << "FROM uniprot_taxonomy "
#else
                            << "FROM ncbi_taxonomy "
#endif
                            << "WHERE id IN " << taxa_clause;
#ifdef USING_POSTGRES
          R = ST.exec(species_name_query.str());
          ST.commit();
#else
          mysqlpp::ResUse taxa_res = species_name_query.use();
#endif
          int taxonId;
          string scientific_name, common_name;
#ifdef USING_POSTGRES
          for (pqxx::result::const_iterator i = R.begin(); i != R.end(); ++i) {
            (*i)[0].to(taxonId);
            if ((*i)[1].is_null()) {
              if (!(*i)[2].is_null()) {
                (*i)[2].to(common_name);
                names_of_taxa[taxonId] = common_name;
              }
            } else {
              (*i)[1].to(scientific_name);
              names_of_taxa[taxonId] = scientific_name;
            }
          }
#else
          while (mysqlpp::Row row = taxa_res.fetch_row()) {
            taxonId = int(row["id"]);
            scientific_name = row["scientific_name"].get_string();
            common_name = row["common_name"].get_string();
            if (scientific_name == "NULL") {
              if (common_name != "NULL") {
                names_of_taxa[taxonId] = common_name;
              }
            } else {
              names_of_taxa[taxonId] = scientific_name;
            }
          }
#endif
          cout << "Getting breadth-first visit order" << endl;
          vector<int> nodes_to_visit;
          getBreadthFirstVisitOrder(*tree, nodes_to_visit);
          cout << "Printing reverse visit order" << endl;
          for (int i = nodes_to_visit.size() - 1; i >= 0; --i) {
            cout << "nodeId: " << nodes_to_visit[i];
            if (tree->hasFather(nodes_to_visit[i])) {
              cout << " ParentId: " << tree->getFatherId(nodes_to_visit[i]);
            }
            if (tree->isLeaf(nodes_to_visit[i])) {
              cout  << " UniProt Accession: "
                    << uniprot_accession_of_sequence_header[
                        retrieved_sequence_header_of_leaf[nodes_to_visit[i]]];
            }
            cout << endl;
          }
          cout << "Getting longest distance in tree" << endl;
            double longest_distance = get_longest_distance_in_tree(*tree,
                                                              nodes_to_visit);
          cout << "Longest distance: " << longest_distance << endl;
          double reciprocal_longest_distance;
            reciprocal_longest_distance = 1.0 / longest_distance;
          if (reciprocal_longest_distance 
              != retrieved_reciprocal_longest_distance) {
            if (retrieved_reciprocal_longest_distance != -1.0) {
              cout << "Retrieved reciprocal_longest_distance " 
                  << retrieved_reciprocal_longest_distance
                  << " is not the reciprocal of " << longest_distance << endl;
            }
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
          find_inparalogs_in_tree(*tree, nodes_to_visit, 
                                  stored_species_of_node);
          cout << "Storing inparalogs in database" << endl;
          if (!putInparalogsInDB(*tree, idMap, db_id_of_left_id,
                                "species", conn, stored_species_of_node)) {
            cout << "Inparalogs were already in database, didn't store" << endl;
          }
          cout << "Retrieving inparalogs from database" << endl;
          getInparalogsFromDB(*tree, nodeIdOfLeftIdMap, tree_id,
                              "species", conn, retrieved_species_of_node);
          cout << "Finding orthologs in tree" << endl;
          map<DBIdentifierT, DBIdentifierT> 
            stored_unique_sequence_header_with_species,
            retrieved_unique_sequence_header_with_species;
          set<DBIdentifierT> species_with_multiple_sequence_headers;
          map<int, map<DBIdentifierT, TreeDistanceInfo<ObjectT,
                                            nullObjectValue> *> *>
            stored_distance_from_sequence_header_with_taxon_to_node,
            retrieved_distance_from_sequence_header_with_taxon_to_node;
          find_orthologs_in_tree(*tree, nodes_to_visit, 
                          stored_sequence_header_of_leaf, 
                          stored_species_of_node,
                          stored_unique_sequence_header_with_species,
                          species_with_multiple_sequence_headers,
                      stored_distance_from_sequence_header_with_taxon_to_node);
          map<DBIdentifierT, DBIdentifierT>::const_iterator
            stored_unique_iter;
          for (stored_unique_iter 
                  = stored_unique_sequence_header_with_species.begin();
              stored_unique_iter != stored_unique_sequence_header_with_species.end();
              ++stored_unique_iter) {
            cout << "Species: " << names_of_taxa[(*stored_unique_iter).first]
                  << " Unique UniProt Accession: " 
                  << uniprot_accession_of_sequence_header[
                                                (*stored_unique_iter).second];
            cout << endl;
          }
          cout << "Finding proximal subtrees" << endl;
          map<int, set<DBIdentifierT> *>
      stored_species_with_multiple_sequence_headers_for_which_node_is_maximal,
    retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal;
          set<int> stored_superorthologous_nodes, 
                  retrieved_superorthologous_nodes;
          map<DBIdentifierT, set<int> *> maximal_nodes_of_sequence_header;
          map<DBIdentifierT, map<DBIdentifierT, set<int> *> *>
              alternative_nearest_sequence_headers;
          find_proximal_subtrees(*tree, nodes_to_visit,
              stored_species_of_node,
              stored_unique_sequence_header_with_species,
              species_with_multiple_sequence_headers,
              stored_distance_from_sequence_header_with_taxon_to_node,
              idMap,
              leaf_of_sequence_header,
              stored_sequence_header_of_leaf,
      stored_species_with_multiple_sequence_headers_for_which_node_is_maximal,
              stored_superorthologous_nodes,
              maximal_nodes_of_sequence_header,
              alternative_nearest_sequence_headers);
          cout << "Printing alternative nearest sequence headers" << endl;
          map<DBIdentifierT, set<int> *>::const_iterator obj_nodes_iter;
          map<DBIdentifierT, map<DBIdentifierT, set<int> *> *>::const_iterator
            obj_objs_nodes_iter;
          set<int>::const_iterator node_iter;
          LeftRightIdsOfNodeIdMap::const_iterator node_left_iter;
          for (obj_objs_nodes_iter 
                = alternative_nearest_sequence_headers.begin();
              obj_objs_nodes_iter
                != alternative_nearest_sequence_headers.end();
              ++obj_objs_nodes_iter) {
            cout << "Sequence header " << (*obj_objs_nodes_iter).first << endl;
            for (obj_nodes_iter = (*obj_objs_nodes_iter).second->begin();
                obj_nodes_iter != (*obj_objs_nodes_iter).second->end();
                ++obj_nodes_iter) {
              cout << "Alt: " << (*obj_nodes_iter).first;
              cout << " LeftIds of Maximal Nodes: ";
              for (node_iter = (*obj_nodes_iter).second->begin();
                  node_iter != (*obj_nodes_iter).second->end();
                  ++node_iter) {
                node_left_iter = idMap.find(*node_iter);
                cout <<  ((*node_left_iter).second).first << " ";
              }
              cout << endl;
            }
          }
          cout << "Putting orthologs in database" << endl;
          if (!putOrthologsInDB(*tree, idMap, db_id_of_left_id,
              "species", "sequence_header", 
              conn,
              stored_unique_sequence_header_with_species,
              stored_distance_from_sequence_header_with_taxon_to_node,
    stored_species_with_multiple_sequence_headers_for_which_node_is_maximal)) {
            cout << "Orthologs were already in database, didn't store" << endl;
          }
          /*
          cout << "Retrieving orthologs from database" << endl;
          getOrthologsFromDB(*tree, nodeIdOfLeftIdMap, tree_id, 
                  "species", "sequence_header",
                  conn,
                  retrieved_unique_sequence_header_with_species,
                  retrieved_distance_from_sequence_header_with_taxon_to_node,
  retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal);
          cout << "Retrieved orthologs from database" << endl;
          */
          map<DBIdentifierT, DBIdentifierT>::const_iterator
            unique_iter;
          for (unique_iter 
                  = stored_unique_sequence_header_with_species.begin();
              unique_iter != stored_unique_sequence_header_with_species.end();
              ++unique_iter) {
            cout << "Species: " << names_of_taxa[(*unique_iter).first]
                  << " Unique UniProt Accession: " 
                  << uniprot_accession_of_sequence_header[
                                                (*unique_iter).second];
            if (retrieved_unique_sequence_header_with_species.find((
                  *unique_iter).first) ==
                retrieved_unique_sequence_header_with_species.end()) {
              cout << " NOT IN DB!!! ";
            } else if (retrieved_unique_sequence_header_with_species[
                        (*unique_iter).first] != (*unique_iter).second) {
              cout << "!=" << (*unique_iter).second;
            }
            cout << endl;
          }
          cout << "Putting " << stored_superorthologous_nodes.size()
              << " superorthologous nodes in database" << endl;
          putSuperorthologousNodesInDB(*tree, idMap, db_id_of_left_id, tree_id,
              "species", "sequence_header", conn, stored_superorthologous_nodes);
          cout << "Retrieving superorthologous nodes from database" << endl;
          getSuperorthologousNodesFromDB(*tree, nodeIdOfLeftIdMap, 
              tree_id, "species", "sequence_header", 
              conn,
              retrieved_superorthologous_nodes);
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
            stored_species_of_node, 
            stored_distance_from_sequence_header_with_taxon_to_node,
            idMap, nodeIdOfLeftIdMap, leaf_of_sequence_header, 
            stored_sequence_header_of_leaf,
            alternative_nearest_sequence_headers,
            duplication_node_of_sequence_headers,
            pair_of_sequence_headers_yielding_duplication_distance,
            pair_of_nearest_sequence_headers_to_maximal_node,
            greatest_distance_of_maximal_descendant);
          cout << "Printing duplication node of sequence headers" << endl;
          map<DBIdentifierT, map<DBIdentifierT, DupInfo *> *>::const_iterator
            seqhdr_seqhdr_dupinfo_iter;
          map<DBIdentifierT, DupInfo *>::const_iterator seqhdr_dupinfo_iter;
          for (seqhdr_seqhdr_dupinfo_iter
                = duplication_node_of_sequence_headers.begin();
              seqhdr_seqhdr_dupinfo_iter 
                != duplication_node_of_sequence_headers.end();
              ++seqhdr_seqhdr_dupinfo_iter) {
            cout << "Sequence Header: " << seqhdr_seqhdr_dupinfo_iter->first
              << endl;
            for (seqhdr_dupinfo_iter 
                = seqhdr_seqhdr_dupinfo_iter->second->begin();
                seqhdr_dupinfo_iter
                != seqhdr_seqhdr_dupinfo_iter->second->end();
                ++seqhdr_dupinfo_iter) {
              cout << "  Sequence Header: " << seqhdr_dupinfo_iter->first;
              node_left_iter = idMap.find(
                        seqhdr_dupinfo_iter->second->getDuplicationNodeId());
              cout << " Leftid of Duplication Node: " 
                    << node_left_iter->second.first << endl;
            }
          }
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
          cout << "Printing tree" << endl;
          printNodeId(*tree, tree->getRootId(), 0, idMap, levelOfNodeIdMap,
              stored_sequence_header_of_leaf,
              retrieved_sequence_header_of_leaf,
              uniprot_accession_of_leaf,
              uniprot_accession_of_sequence_header,
              names_of_taxa,
              stored_species_of_node,
              retrieved_species_of_node,
              stored_distance_from_sequence_header_with_taxon_to_node,
/*
              retrieved_distance_from_sequence_header_with_taxon_to_node,
*/
      stored_species_with_multiple_sequence_headers_for_which_node_is_maximal,
/*
    retrieved_species_with_multiple_sequence_headers_for_which_node_is_maximal,
*/
              stored_superorthologous_nodes,
              retrieved_superorthologous_nodes,
              duplication_node_of_sequence_headers,
              pair_of_sequence_headers_yielding_duplication_distance,
              pair_of_nearest_sequence_headers_to_maximal_node,
              greatest_distance_of_maximal_descendant);
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
          for (obj_nodes_iter = maximal_nodes_of_sequence_header.begin();
              obj_nodes_iter != maximal_nodes_of_sequence_header.end();
              ++obj_nodes_iter) {
            delete (*obj_nodes_iter).second;
          }
          for (obj_nodes_iter = maximal_nodes_of_sequence_header.begin();
              obj_nodes_iter != maximal_nodes_of_sequence_header.end();
              ++obj_nodes_iter) {
            delete (*obj_nodes_iter).second;
          }
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
          for (dist_attr_node_iter = 
              stored_distance_from_sequence_header_with_taxon_to_node.begin();
                dist_attr_node_iter != 
              stored_distance_from_sequence_header_with_taxon_to_node.end();
                ++dist_attr_node_iter) {
            for (attr_node_iter = ((*dist_attr_node_iter).second)->begin();
                  attr_node_iter != ((*dist_attr_node_iter).second)->end();
                  ++attr_node_iter) {
              delete (*attr_node_iter).second;
            }
            delete (*dist_attr_node_iter).second;
          }
          for (dist_attr_node_iter = 
            retrieved_distance_from_sequence_header_with_taxon_to_node.begin();
                dist_attr_node_iter != 
            retrieved_distance_from_sequence_header_with_taxon_to_node.end();
                ++dist_attr_node_iter) {
            for (attr_node_iter = ((*dist_attr_node_iter).second)->begin();
                  attr_node_iter != ((*dist_attr_node_iter).second)->end();
                  ++attr_node_iter) {
              delete (*attr_node_iter).second;
            }
            delete (*dist_attr_node_iter).second;
          }
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
