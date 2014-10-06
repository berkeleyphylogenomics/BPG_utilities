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

#include "AlignmentStatistics.h"
#include "TreeSVG.h"
#include "AlignmentSVG.h"

#define DEBUG true

// Methods
Fasta getFasta(Fasta *fa, stack<string> names);
float averageBootstrap(Node *n);
bool cutNode(Tree *tree, Node *node, Fasta *fasta, int minbs, int minpid);
void computePercentId(Tree *tree, Fasta *msa);
void saveSubfam(Tree *tree, Fasta *msa, Node *node, int subfamCount);
void debug(Node *node, int *nodeCount);
void usage();

// Globals (!)
int subfamCount = -1;
string outtree = "";
ofstream summaryofs;
ofstream treesofs;
ofstream msaofs;
ofstream treeofs;

/*
 * kerp.cpp
 *
 * Usage:  kerp <minbs> <minpid> <tree> <msa>
 *   minbs:  The minimum bootstrap value as an integer. The maximum is the
 *           number of datasets used when bootstrapping
 *   minpid: The minimum percent ID as an integer
 *   tree:   The tree file in Newick format (must have bootstrap values)
 *   msa:    The alignment in FASTA or A2M format
 *
 * Output: Four files:
 *   kerf.summary: A summary of subfamilies
 *   kerf.seqs:    The sequences in each subfamily, separated by '#'
 *   kerf.trees:   The subfamily trees, separated by '#'
 *   kerf.tre:     A multifurcating Newick tree containing all subfamilies
 *
 * Author:  John H. Lee <johnhlee@berkeley.edu>
 *          Berkeley Phylogenomics Group -- http://phylogenomics.berkeley.edu
 */ 
int main(int argc, char *argv[])
{
    // Check args
    if (argc < 5)
    {
        usage();
        return(0);
    }

    // Read input params
    int minbs = atoi(argv[1]);
    cout << "Minimum bootstrap value: " << minbs << endl;
    int minpid = atoi(argv[2]);
    cout << "Minimum percent id: " << minpid << endl;
    string treefn = argv[3];
    string msafn  = argv[4];

    cout << "Reading files...";
    cout.flush();

    // Read the tree
    Tree tree;
    tree.readFile(treefn);

    // Read the alignment
    Fasta msa;
    msa.readFile(msafn);

    // Open and configure output streams
    summaryofs.open("kerf.summary");
    summaryofs.precision(2);
    summaryofs.setf(ios::right);
    summaryofs.setf(ios::fixed);
    treesofs.open("kerf.trees");
    treeofs.open("kerf.tre");
    msaofs.open("kerf.seqs");

    cout << "done" << endl;
    cout.flush();

//    cout << tree.toNewickString() << endl;

    cout << "Computing percent ID...";
    cout.flush();

    // Compute percent ID at each node
    computePercentId(&tree, &msa);

    cout << "done" << endl;
    cout.flush();

    // Build the final multifurcating tree here
    outtree = "(";

    // Print the summary header.  The rest is written in saveSubfam()
    summaryofs << "Subfamily    Sequences    Percent Identity    Bootstrap";
    summaryofs << endl;
    summaryofs << "---------    ---------    ----------------    ---------";
    summaryofs << endl;

    cout << "Identifying subfamilies...";
    cout.flush();

    // Cut the tree
    cutNode(&tree, tree.getRoot(), &msa, minbs, minpid);

    cout << "done" << endl;
    cout.flush();

    // The last char of outtree is a comma.  Replace with ) and close with ;
    outtree[outtree.size() - 1] = ')';
    outtree += ";";

    // Write the tree file
    treeofs << outtree.c_str() << endl;

    // Close files
    summaryofs.close();
    treesofs.close();
    treeofs.close();
    msaofs.close();

    cout << "Writing results...";
    cout.flush();

    // Write the full SVG tree (soon, with highlighted subfams)
    //Tree restree;
    //restree.readFile(treefn.c_str());
    TreeSVG svg;
    ofstream svgofs;
    svgofs.open("kerf.svg");
    //svgofs << svg.toString(&restree) << endl;
    svgofs << svg.toString(&tree) << endl;
    svgofs.close();
/*
    AlignmentSVG msasvg;
    ofstream msasvgofs;
    msasvgofs.open("kerf_msa.svg");
    msasvgofs << msasvg.toString(&msa) << endl;
    msasvgofs;

    system("svg2png kerf_msa.svg >/dev/null 2>&1");
*/
    system("svg2png kerf.svg >/dev/null 2>&1");

    cout << "done" << endl;

    // Fin
    return(0);
}

/*
 * Traverse the tree in inorder, computing percent identity at each node.
 * The result is stored in each Node.
 */
void computePercentId(Tree *tree, Fasta *msa)
{
    stack<Node*> s;

    Node *n = tree->getRoot();

    int nseqs = msa->count();

    float *percentIdTable = new float[nseqs*nseqs];
    for (int i = 0; i < nseqs*nseqs; i++)
    {
        percentIdTable[i] =-1;
    }
    
    AlignmentStatistics astats;
    astats.allowGaps(true);
    
    for ( ; ; )
    {
        if (n != NULL)
        {
            s.push(n);
            n = n->getChild(0);
        }
        else
        {
            if (s.size() == 0)
            {
                break;
            }

            n = s.top();
            s.pop();

            // Begin visit
            cout << "Visiting node";
            cout << " name " << n->getName();
            cout << " childCount " << n->getChildCount();
            cout << " Level " << n->getLevel();
            cout << " Index " << n->getIndex() << endl;

            if (n->isInternal())
            {
                float minID = 100;
                //stack<string> names = tree->leafNames(n);
                vector<Node*> leaves = tree->leaves(n);
                for (int i = 0; i < leaves.size(); i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        Node *node_i = leaves.at(i);
                        Node *node_j = leaves.at(j);
                        int index_i = node_i->getIndex();
                        int index_j = node_j->getIndex();
                        
                        // if we haven't calculated the %ID for this pair yet
                        if (percentIdTable[index_j * nseqs + index_i ] == -1)
                        {
                            Sequence seq_i = msa->getSequence(node_i->getName());
                            Sequence seq_j = msa->getSequence(node_j->getName());
                            if (seq_i.getId() == "" || seq_j.getId() == "")
                            {
                                cerr << "There is a mismatch between the IDs in the tree and the alignment!" << endl;
                            }
                            float percentId = astats.pairPercentId(&seq_i, &seq_j);
                            percentIdTable[index_i * nseqs + index_j] = percentId;
                            percentIdTable[index_j * nseqs + index_i] = percentId;
                        }

                        if (percentIdTable[index_i * nseqs + index_j] < minID)
                        {
                            minID = percentIdTable[index_i * nseqs + index_j];
                        }
                    }
                }
                //Fasta nodefa = getFasta(msa, names);
                //AlignmentStatistics astats;
                //astats.allowGaps(false); //use new method for %id calculation
                //astats.compute(nodefa);
                n->setPercentId(minID);
                cout << "Set minimum percent id " << minID;
                cout << " over " << leaves.size() << " leaves ";
                cout << " from " << n->getChildCount() << " children." << endl;
            }

            // End visit

            //cout << "Just visited node " << n->getName() << ":" << n->getIndex() << endl;
            
            n = n->nextSibling();
            /*
            Node *parent = n->getParent();
            if (parent != NULL)
            {
                int childIndex = parent->getChildIndex(n);
                n = parent->getChild(childIndex + 1);
                //cout << "Now visiting node " << n->getIndex() << "." << endl;
            }
            else
            {
                n = NULL;
            }*/
        }
    }
    delete [] percentIdTable;
} 

/*
 * Recursively attempt to cut the tree from node.  If the node fails
 * our PID and bootstrap tests, recurse to children.  Otherwise cut and
 * save the subfam.
 */
bool cutNode(Tree *tree, Node *node, Fasta *msa, int minbs, int minpid)
{
    if (node == NULL)
    {
        cout << "node is NULL!" << endl;
        return false;
    }

    if (node->getSubfam() != -1)
    {
        cout << "subfam already set" << endl;
        return false;
    }

    // If this is an external node, we have no choice but to make it a subfam
    // DB: with bottum up, we return true if the lower node passes the bs & pctid tests
    if (node->isExternal())
    {
//        saveSubfam(tree, msa, node, ++subfamCount);
        return true;
    }

/* Bottom-up */
    int children = node->getChildCount();
    bool *test = new bool[children];
    bool allGood = true;
    
    for (int i = 0; i < children; i++)
    {
        test[i] = cutNode(tree, node->getChild(i), msa, minbs, minpid);
        if (! test[i]) allGood = false;
    }
//    bool testRight = cutNode(tree, node->getRight(), msa, minbs, minpid);
//    bool testLeft = cutNode(tree, node->getLeft(), msa, minbs, minpid);

    // if both left and right nodes pass the test, we can try to merge them into a subfam.
    if (allGood)
    {
        if (node->getPercentId() >= minpid || 
            (minbs > 0 && node->getBootstrapInt() >= minbs 
                       && node->getParent() != NULL ))
        {
            delete test;
            return true;
        }
        else
        {
            for (int i = 0; i < children; i++)
            {
                saveSubfam(tree, msa, node->getChild(i), ++subfamCount);
            }
            // cut here and define left and right as two subfamilies
            //saveSubfam(tree, msa, node->getRight(), ++subfamCount);
            //saveSubfam(tree, msa, node->getLeft(), ++subfamCount);
            delete test;
            return false;
        }
    }
    else
    {
        for (int i = 0; i < children; i++)
        {
            if (test[i])
            {
                saveSubfam(tree,msa, node->getChild(i), ++subfamCount);
            }
        }

        delete test;
        return false;
    }
    /*
        if ( testRight )
        {
            saveSubfam(tree, msa, node->getRight(), ++subfamCount);
        }
        else if ( testLeft )
        {
            saveSubfam(tree, msa, node->getLeft(), ++subfamCount);
        }

        return false;
        
    } */

/*
    if (cutNode(tree, node->getRight(), msa, minbs, minpid) ||
        cutNode(tree, node->getLeft(), msa, minbs, minpid))
    {
        node = node->getParent();

        if (node == NULL)
        {
            cout << "Moving this node to parent" << endl;
            return false;
        }
    }

    if (node->getPercentId() > minpid || node->getBootstrapInt() > minbs)
    {
        saveSubfam(tree, msa, node, ++subfamCount);
        return true;
    }
*/

    // If the node has low PID and BS, recurse to the children
/*  Top-down  */
/*    if (node->getPercentId() < minpid || node->getBootstrapInt() < minbs)
    {
        cutNode(tree, node->getRight(), msa, minbs, minpid);
        cutNode(tree, node->getLeft(), msa, minbs, minpid);
    }
    else    // Else, save the subfam
    {
        saveSubfam(tree, msa, node, ++subfamCount);
        return true;
    }
*/
    return false;
}

/*
 * Write a subfamily to our output files and reposition the node's sibling
 */
void saveSubfam(Tree *tree, Fasta *msa, Node *node, int subfamCount)
{
    // Get Fasta for all leafs in this node
    stack<string> names = tree->leafNames(node);
    Fasta fa = getFasta(msa, names);

    // Get Newick for the node
    string newick = tree->newickStringForNode(node);

    // Add this tree to the main, multifurcating output tree
    ostringstream nodeinfo;
    nodeinfo << node->getParent()->getPercentId() << "_" <<
        node->getParent()->getBootstrap() << ":" <<
        node->getParent()->getBootstrap();

    outtree += "(" + newick + ")" + nodeinfo.str() + ",";

    // Write the summary
    summaryofs << setw(9) << subfamCount;
    summaryofs << setw(13) << fa.count();
    summaryofs << setw(20) << node->getPercentId();
    summaryofs << setw(13) << node->getBootstrap();
    summaryofs << endl;

    // Make sure we have valid Newick
    if (node->isExternal())
    {
        newick = "(" + newick + ");"; 
    }
    else
    {
        newick += ";";
    }

    // Write the tree and msa to our output files
    treesofs << "\%subfamily " << subfamCount << endl;
    treesofs << newick;
    treesofs << endl;

    msaofs << "\%subfamily " << subfamCount << endl;
    msaofs << fa.toFasta();
    msaofs << endl;

    //XXX If we want to write out separate files for each subfam
    if (DEBUG)
    {
        stringstream treefn;
        treefn << "subfam" << subfamCount << ".tre";
        stringstream msafn;
        msafn << "subfam" << subfamCount << ".fa";

        ofstream fos;
        string data = newick;
        fos.open(treefn.str().c_str());
        fos.write(data.c_str(), data.size());
        fos.close();

        data = fa.toFasta();
        fos.open(msafn.str().c_str());
        fos.write(data.c_str(), data.size());
        fos.close();
    }

    node->setSubfam(subfamCount);
/*
    // Reassign our sibling to our grandparent
    Node *sibling = node->getSibling();
    if (sibling == NULL) { return; }

    Node *parent  = node->getParent();
    if (parent == NULL) { return; }

    Node *grandParent = parent->getParent();

    if (grandParent == NULL)
    {
        tree->setRoot(parent);
    }
    else
    {
        sibling->setParent(grandParent);

        if (grandParent->getLeft() == parent)
        {
            grandParent->setLeft(sibling);
        }
        else if (grandParent != NULL)
        {
            grandParent->setRight(sibling);
        }
    }
*/
}

/*
 * Get Fasta for a list of names (analogous to fastacmd)
 */
Fasta getFasta(Fasta *fa, stack<string> names)
{
    Fasta subfam;

    while (! names.empty())
    {
        string name = names.top();
        names.pop();

        Sequence s = fa->getSequence(name);

        if (s.getId() == "")
        {
            cerr << "ERROR: null sequence for name " << name << endl;
        }
        else
        {
            subfam.addSequence(fa->getSequence(name));
        }
    }

    return subfam;
}

/*
 * Compute average bootstrap for a node.  XXX Unused
 */
float averageBootstrap(Node *n)
{
    // Traverse the tree in inorder

    int nodeCount    = 0;
    int bootstrapSum = 0;

    stack<Node*> s;

    for ( ; ; )
    {
        if (n != NULL)
        {
            s.push(n);
            n = n->getChild(0);
        }
        else
        {
            if (s.size() == 0)
            {
                break;
            }

            n = s.top();
            s.pop();

            // Begin visit

            nodeCount++;
            bootstrapSum += n->getBootstrapInt();

            // End visit

            Node *parent = n->getParent();
            if (parent != NULL)
            {
                int childIndex = parent->getChildIndex(n);
                n = n->getChild(childIndex + 1);
            }
//            n = n->getRight();
        }
    }

    return (float) bootstrapSum / nodeCount;
}

/*
 * Debug
 */
void debug(Node *node, int *nodeCount)
{
    if (node != NULL)
    {
        for (int i = 0; i < node->getChildCount(); i++)
        {
            debug(node->getChild(i), nodeCount);
        }
        /*
        if (node->hasRight())
        {
            debug(node->getRight(), nodeCount);
        }

        if (node->hasLeft())
        {
            debug(node->getLeft(), nodeCount);
        }

        */
        cout << "Node:     " << *nodeCount << endl;
        cout << "    Name: " << node->getName() << endl;
        cout << "    PID:  " << node->getPercentId() << endl;
        cout << "    BS:   " << node->getBootstrap() << endl;

        for (int i = 0; i < node->getChildCount(); i++)
        {
            cout << "    " << i << ":    " << node->getChild(i)->getName() << endl;
        }
        /*
        if (node->hasRight())
        {
            cout << "    R:    " << node->getRight()->getName() << endl;
        }

        if (node->hasLeft())
        {
            cout << "    L:    " << node->getLeft()->getName() << endl;
        }

        */
        *nodeCount = *nodeCount + 1;
    }
}

/*
 * Print usage info
 */
void usage()
{
    cout << endl;
    cout << "kerf 1.0 -- Cut a Newick tree on bootstrap and percent identity";
    cout << endl << endl;
    cout << "Usage: kerf <minbs> <minpid> <treefn> <msafn>" << endl;
    cout << endl;
    cout << "Parameters:" << endl;
    cout << "  <minbs>:   Minimum bootstrap value as an integer" << endl;
    cout << "  <minpid :  Minimum percent identity as an integer" << endl;
    cout << "  <treefn>:  Newick tree file with bootstrap values" << endl;
    cout << "  <msafn>:   Alignment in FASTA or A2M format" << endl;
    cout << endl;
    cout << "Berkeley Phylogenomics Group -- http://phylogenomics.berkeley.edu";
    cout << endl;
    cout << endl;
}

