/* connection to postgres db */
#include <pqxx/pqxx>

/* user-input options parsing */
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include <exception>
using std::exception;

#include <map>
using std::map;

#include <vector>
using std::vector;

#include <Phyl/Newick.h>
using bpp::Newick;

#include <Phyl/Tree.h>
using bpp::Tree;

#include "modified_preorder_tree_traversal.h"
#include "put_left_right_ids_in_db.h"
#include "get_taxon_of_leaf_map.h"
#include "get_breadth_first_visit_order.h"
#include "get_longest_distance_in_tree.h"
#include "find_inparalogs_in_tree.h"
#include "put_inparalogs_in_db.h"
#include "find_orthologs_in_tree.h"
#include "find_proximal_subtrees.h"
#include "put_orthologs_in_db.h"
#include "put_superorthologous_nodes_in_db.h"


/*
main routine for the find_ortholog program. The program takes two options:

1. --book book_id, where book_id is a bpg_accession number (just the number, with zeros okay) currently in the database
2. --tree method, where method is nj, ml, or sciphy. This tree must have been created and must be in the appropriate place
                                                     in the directory hierarchy

The program connects to the database, looks up the relevant information about the book, finds orthologs
from the given tree, and finally updates the relevant database tables.
*/
int main(int argc, char* argv[])
{
    /* DB access parameters */
    string host = "db";
    string database = "pfacts003_test";
    string user = "webuser";
    string password = "w3zx19Ko";

    string dbaccess = "host=" + host + " dbname=" + database + " user=" + user + " password=" + password;

    /* location of flat files for phylofacts */
    string filedir = "/clusterfs/ohana/bpg/phylofacts";

    if(database.compare("pfacts003_test") == 0)
        filedir = "/clusterfs/ohana/bpg/phylofacts_test";

    /* parse user options */
    try 
    {
        po::options_description desc("Allowed options");
    
        desc.add_options()
            ("help", "produce help message")
            ("book", po::value<string>(), "book within which to find orthologs")
            ("tree", po::value<string>(), "tree method to find orthologs for (nj, ml, sciphy), default is nj")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        /* print help message */
        if(vm.count("help"))
        {
            cout << desc << endl;
            return 1;
        }

        /* read tree method */
        string tree_method = "nj";

        if(vm.count("tree"))
            tree_method = vm["tree"].as<string>();

        if(tree_method.compare("nj") != 0 and tree_method.compare("ml") != 0 and tree_method.compare("sciphy") != 0)
        {
            cout << "tree must be one of: nj, ml, sciphy" << endl;
            return 1;
        }

        /* read book id */
        if(! vm.count("book"))
        {
            cout << "book id required to run this program" << endl;
            return 1;
        }

        string book_id = vm["book"].as<string>();

        /* done parsing user input, run the thing, already */

        /* get db connection */
        pqxx::connection conn(dbaccess);
        pqxx::transaction<> xaction(conn);

        cout << "checking that book id is valid...";

        /* make sure family (book id) is in db  */
        pqxx::result res = xaction.exec("SELECT id FROM family where id=" + book_id);

        if(res.empty())
        {
            cout << "Book " << book_id << " not found in the database." << endl;
            return 1;
        }

        cout << "looks good." << endl;
        cout << "checking that requested tree is there...";

        /* make sure requested tree is in db */
        res = xaction.exec("SELECT tree FROM family_tree WHERE family_id=" + book_id + " AND tree_method='" + tree_method + "'");

        if(res.empty())
        {
            cout << "Tree " << tree_method << " not found for book " << book_id << endl;
            return 1;
        }

        Tree* tree = 0;

        istringstream treestream(res[0]["tree"].c_str());

        Newick* reader = new Newick(false);
        tree = reader->read(treestream);
        delete reader;

        if(tree == 0)
        {
            cout << "Error reading tree" << endl;
            return 1;
        }

        cout << "done." << endl;

        /* Now we have the tree, get its left and right ids */
        int max_right_id;
        LeftRightIdsOfNodeIdMap idMap;
        map<int,int> nodeIdOfLeftIdMap;
        map<int,int> levelOfNodeIdMap;
 
        cout << "mapping node ids to left and right ids...";

        max_right_id = findLeftRightIds(*tree, idMap, nodeIdOfLeftIdMap, levelOfNodeIdMap, tree->getRootId(), 1, 0);

        cout << "done." << endl;

        /* put tree into db */
        cout << "putting tree into database...";

        xaction.exec("INSERT INTO tree (family_id, method, is_rsd_rooted) VALUES (" + book_id + ", \'" + tree_method + "\', FALSE)");

        /* get tree id */
        res = xaction.exec("SELECT last_value FROM tree_id_seq");

        string tree_id = res[0][0].c_str();

        cout << "done." << endl;

        /* put left and right ids into database  */
        cout << "putting left and right ids into database...";

        map<int,int> protein_sequence_of_leaf;
        map<int,int> db_id_of_left_id;

        putLeftRightIdsInDB(*tree, idMap, nodeIdOfLeftIdMap, levelOfNodeIdMap, atoi(tree_id.c_str()), book_id, db_id_of_left_id, protein_sequence_of_leaf, xaction);

        cout << "done." << endl;

        /* get species identifiers for all leaves */
        map<int,int> species_of_node;

        cout << "getting species for all leaves...";

        getTaxonOfLeafMap(protein_sequence_of_leaf, species_of_node, xaction);

        cout << "done." << endl;   

        /* get breadth-first order */
        vector<int> nodes_to_visit;

        cout << "getting breadth-first visit order...";

        getBreadthFirstVisitOrder(*tree, nodes_to_visit);
    
        cout << "done." << endl;

        /* get reciprocal longest distance */
        cout << "getting reciprocal of longest distance in tree...";

        double longest_distance = get_longest_distance_in_tree(*tree, nodes_to_visit);
        double reciprocal_longest_distance = 1.0 / longest_distance;

        ostringstream query;

        query << "UPDATE tree SET reciprocal_longest_distance=" << reciprocal_longest_distance
              << " WHERE id=" << tree_id;

        xaction.exec(query.str());

        cout << "done." << endl;

        /* find inparalogs in the tree */
        
        cout << "finding inparalogs in tree...";

        find_inparalogs_in_tree(*tree, nodes_to_visit, species_of_node);

        cout << "done." << endl;

        /* store inparalogs in database */

        cout << "putting inparalogs in database...";

        putInparalogsInDB(*tree, idMap, db_id_of_left_id, "species", species_of_node, xaction);

        cout << "done." << endl;

        /* find orthologs in tree */
        
        cout << "finding orthologs in tree...";

        map<int,int> unique_protein_sequence_with_species;
        set<int> species_with_multiple_protein_sequences;
        map<int, map<int, TreeDistanceInfo<ObjectT, nullObjectValue> *> *> distance_from_protein_sequence_with_taxon_to_node;

        find_orthologs_in_tree(*tree, nodes_to_visit, protein_sequence_of_leaf, species_of_node,
                               unique_protein_sequence_with_species, species_with_multiple_protein_sequences,
                               distance_from_protein_sequence_with_taxon_to_node);

        cout << "done." << endl;

        /* fing proximal subtrees */

        cout << "finding proximal subtrees...";

        map<int, set<int> *> species_with_multiple_protein_sequences_for_which_node_is_maximal;
        set<int> superorthologous_nodes;
        map<int, set<int> *> maximal_nodes_of_protein_sequence;

        find_proximal_subtrees(*tree, nodes_to_visit, species_of_node,
                               unique_protein_sequence_with_species,
                               species_with_multiple_protein_sequences,
                               distance_from_protein_sequence_with_taxon_to_node,
                               species_with_multiple_protein_sequences_for_which_node_is_maximal,
                               superorthologous_nodes, maximal_nodes_of_protein_sequence);

        cout << endl << "  there are " << unique_protein_sequence_with_species.size()
                     << " species with unique protein sequences in the tree." << endl;

        cout << "  there are " << species_with_multiple_protein_sequences_for_which_node_is_maximal.size()
             << " species with multiple protein sequences in the tree." << endl;

        cout << "done." << endl;

        /* put orthologs in db */

        cout << "putting orthologs in database...";

        putOrthologsInDB(*tree, idMap, db_id_of_left_id,
                         "species", "protein_sequence",
                         unique_protein_sequence_with_species,
                         distance_from_protein_sequence_with_taxon_to_node,
                         species_with_multiple_protein_sequences_for_which_node_is_maximal,
                         xaction);

        cout << "done." << endl;

        /* put superorthologous nodes in db */

        cout << "putting " << superorthologous_nodes.size()
             << " superorthologous nodes in database...";

        putSuperorthologousNodesInDB(*tree, idMap, db_id_of_left_id, 
                                     atoi(tree_id.c_str()),
                                     "species", "protein_sequence",
                                      superorthologous_nodes, xaction);

        cout << "done." << endl;

        /* clean up by freeing memory */

        map<int, set<int> *>::const_iterator nodes_obj_iter;
          
        for (nodes_obj_iter = maximal_nodes_of_protein_sequence.begin();
              nodes_obj_iter != maximal_nodes_of_protein_sequence.end();
              ++nodes_obj_iter) 
        {
            delete (*nodes_obj_iter).second;
        }
          
        map<int, map<int, TreeDistanceInfo<ObjectT, nullObjectValue> *> *>::const_iterator dist_attr_node_iter;
        map<int, TreeDistanceInfo<ObjectT, nullObjectValue> *>::const_iterator attr_node_iter;
          
        for (dist_attr_node_iter
                  = distance_from_protein_sequence_with_taxon_to_node.begin();
                dist_attr_node_iter
                  != distance_from_protein_sequence_with_taxon_to_node.end();
                ++dist_attr_node_iter) 
        {
            for (attr_node_iter = ((*dist_attr_node_iter).second)->begin();
                  attr_node_iter != ((*dist_attr_node_iter).second)->end();
                  ++attr_node_iter) 
            {
                delete (*attr_node_iter).second;
            }
        
            delete (*dist_attr_node_iter).second;
        }
          
        delete tree;
 
        /* finally, commit the transaction */
        xaction.commit();

        cout << "database transaction committed. I'm done." << endl;
    }
    catch(exception& e)
    {
        cout << "error: " << e.what() << endl;
        return 1;
    }

    return 0;
}

