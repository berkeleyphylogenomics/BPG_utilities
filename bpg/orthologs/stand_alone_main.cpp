#include <boost/program_options.hpp>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "modified_preorder_tree_traversal.h"
#include "fill_sequence_header_leaf_maps.h"
#include "get_taxon_of_leaf_map_from_file.h"
#include "get_breadth_first_visit_order.h"
#include "find_inparalogs_in_tree.h"
#include "find_orthologs_in_tree.h"
#include "find_proximal_subtrees.h"
#include "find_duplication_distances.h"
#include "write_phogs_at_threshold.h"
#include "write_phyloxml_file.h"

int main(int ac, char* av[])
{
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("tree", po::value<string>(), "Newick format tree file")
      ("taxa", po::value<string>(), "File with gene,taxon pairs")
      ("xml", po::value<string>(), "File for complete output in PhyloXML format")
      ("threshold", po::value<double>(), "Threshold at which to print PHOGs")
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    if (vm.count("tree") && vm.count("taxa")) {
          string tree_file = vm["tree"].as<string>();
          string taxa_file = vm["taxa"].as<string>();
          bpp::Newick * newickReader = new bpp::Newick(false);
          bpp::Tree * tree;
          cout << "Reading tree from " << tree_file << endl;
          tree = newickReader->read(tree_file);
          int max_right_id;
          LeftRightIdsOfNodeIdMap idMap;
          map<int, int> nodeIdOfLeftIdMap;
          map<int, int> levelOfNodeIdMap;
          cout << "Doing modified preorder tree traversal" << endl;
          max_right_id = findLeftRightIds(*tree, idMap, nodeIdOfLeftIdMap,
                                          levelOfNodeIdMap,
                                          tree->getRootId(), 1, 0);
          map<string, DBIdentifierT> sequence_header_of_leaf_name;
          vector<string> leaf_names;
          map<int, DBIdentifierT> sequence_header_of_leaf;
          map<DBIdentifierT, int> leaf_of_sequence_header;
          cout << "Creating leaf name mapping " << endl;
          fillSequenceHeaderLeafMaps(*tree, 
                                    leaf_names,
                                    sequence_header_of_leaf_name,
                                    sequence_header_of_leaf,
                                    leaf_of_sequence_header);
          map<int, DBIdentifierT> species_of_node;
          vector<string> taxon_names;
          map<string, DBIdentifierT> taxon_of_taxon_name;
          cout << "Retrieving species at leaves " << endl;
          getTaxonOfLeafMapFromFile(taxa_file, leaf_of_sequence_header,
                            sequence_header_of_leaf_name, species_of_node,
                            taxon_names, taxon_of_taxon_name);
          vector<int> nodes_to_visit;
          cout << "Constructing breadth-first visit order " << endl;
          getBreadthFirstVisitOrder(*tree, nodes_to_visit);
          cout << "Finding inparalogs in tree" << endl;
          find_inparalogs_in_tree(*tree, nodes_to_visit, species_of_node);
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
              << " species with unique genes in the tree." << endl;
          cout << "There are " 
    << species_with_multiple_sequence_headers_for_which_node_is_maximal.size()
            << " nodes which are maximal for species with multiple genes"
            << " in the tree." << endl;
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
          cout << "Done computing orthologs" << endl;
          if (vm.count("threshold")) {
            double threshold = vm["threshold"].as<double>();
            writePHOGsAtThreshold(*tree, idMap,
                          duplication_node_of_sequence_headers,
                          pair_of_sequence_headers_yielding_duplication_distance,
                          pair_of_nearest_sequence_headers_to_maximal_node,
                          greatest_distance_of_maximal_descendant,
                          threshold);
          }
          if (vm.count("xml")) {
            string xmlfile = vm["xml"].as<string>();
            write_PhyloXML_File(*tree, idMap, leaf_of_sequence_header,
                          duplication_node_of_sequence_headers,
                          pair_of_sequence_headers_yielding_duplication_distance,
                          pair_of_nearest_sequence_headers_to_maximal_node,
                          greatest_distance_of_maximal_descendant,
                          species_of_node,
                          taxon_names,
                          xmlfile);
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
    } else {
      cout << "Both --tree and --taxa need to be specified." << endl;
      cout << "The --tree option should specify the path to a "
          << "Newick format file containing a gene tree."
          << endl;
      cout << "The --taxa option should specify the specify the path to a "
           << "file specifying the taxon of each gene in the gene tree."
            << endl;
      cout << "Each line of the taxa file should have a gene identifier, "
          << "a tab character, and a taxon identifier, like this:"
          << endl << endl;
      cout << "P0ABB0\t83333" << endl << endl;
      cout << "or this:" << endl << endl;
      cout << "ATPA_ECOLI\t83333" << endl << endl;
      cout << "or this:" << endl << endl;
      cout << "AP_004053.1\t\"Escherichia coli (strain K12)\" " << endl << endl;
      cout << "The gene identifiers can have any format, but they must match "
          << "exactly the entire text of the ones in the gene tree." << endl;
      cout << "The taxon identifiers can also have any format." << endl;
      cout << "If bacterial species and strains are both present, "
          << "that is okay, as long as each is uniquely identified." << endl;
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
