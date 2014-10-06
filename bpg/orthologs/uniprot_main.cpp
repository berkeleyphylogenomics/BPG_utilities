#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "find_orthologs_in_uniprot_tree.h"

int main(int ac, char* av[])
{
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("tree", po::value<string>(), "tree within which to find orthologs")
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    if (vm.count("tree")) {
      string tree_path = vm["tree"].as<string>();
      bpp::Newick * newickReader = new bpp::Newick(false);
      bpp::Tree * tree;
      try {
        tree = newickReader->read(tree_path);
        find_orthologs_in_tree(tree);
      } catch (exception &e) {
        cout << "Error when reading tree." << endl;
      }
      delete tree;
      delete newickReader;

    } else {
      cout << "Tree file was not set.\n";
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
