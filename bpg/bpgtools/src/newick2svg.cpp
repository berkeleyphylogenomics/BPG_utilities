/****************************************************************************
 * Copyright (c) 1998-2005.
 * The Regents of the University of California (Regents).  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this  software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above  copyright notice, this paragraph and the following
 * two paragraphs appear in all copies, modifications, and distributions.  For
 * commercial licensing opportunities, contact:
 *
 *     The Office of Technology Licensing, UC Berkeley
 *     2150 Shattuck Avenue, Suite 510
 *     Berkeley, CA 94720-1620
 *     USA
 *     +1 510 643 7201
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,  INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING  LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS  DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 ****************************************************************************/

#include <iomanip>
#include <stack>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

#include "TreeSVG.h"

#define DEBUG true

// Methods
void usage();

// Globals (!)

/*
 * newick2svg.cpp
 *
 * Usage:  newick2svg <tree>
 *   tree:   The tree file in Newick format
 *
 * Output: SVG on stdout
 *
 * Author:  John H. Lee <johnhlee@berkeley.edu>
 *          Berkeley Phylogenomics Group -- http://phylogenomics.berkeley.edu
 */ 
int main(int argc, char *argv[])
{
  // Check args
  if (argc < 2) //argc always >= 1 since program is first argument
    {
      usage();
      return(0);
    }
  string subfam = "";
  string treefn = argv[1];
  string imagemap = "";
  bool subfamMode = false;
  bool imageMode = false;
  int imageWidth = 320;
    
    /*if (argv[1] == "-s")
    {
        subfam = argv[2];
	cout << subfam << endl;
	cout << argv[3] << endl;
        treefn = argv[3];
	subfamMode = true;
    }*/
  if (argc >= 3)
    {
      subfam = argv[2];
      subfamMode = true;
    }
  if(argc >= 4) 
    {
      imagemap = argv[3];
      imageMode = true;
    }
    if(argc == 5)
      {
	imageWidth = atoi(argv[4]);
      }
    // Read the tree
    Tree tree;
    tree.readFile(treefn);
    tree.pruneSingles(tree.getRoot());
    if(subfamMode) 
      {
	cout << "Adding subfamilies to tree\n";
	tree.readSubfams(subfam);
      }
    TreeSVG svg(imageWidth, true);

    string svgfn = treefn + ".svg";

    cout << "Converting " << treefn << " to " << svgfn << "...";
    ofstream svgofs;
    svgofs.open(svgfn.c_str());
    svgofs << svg.toString(&tree) << endl;
    svgofs.close();
    if(imageMode)
      {
	cout << endl << "Creating imagemap...";	
	svg.printImageMap(&tree,imagemap);
      }
    cout << endl;

    return(0);
}

void usage()
{
    cout << endl;
    cout << "newick2svg 1.0 -- Convert Newick to SVG" << endl;
    cout << endl;
    cout << "Usage: newick2svg <treefn> [subfamfn]" << endl;
    cout << endl;
    cout << "Parameters:" << endl;
    cout << "  <treefn>: Newick tree file name" << endl;
    cout << "  [subfamfn]: subfamily file name, for tree coloring (optional)" << endl;
    cout << endl;
    cout << "Berkeley Phylogenomics Group -- http://phylogenomics.berkeley.edu";
    cout << endl;
    cout << endl;
}

