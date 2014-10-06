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

#ifndef TREE_H
#define TREE_H

#include <string>
#include <math.h>
#include <stack>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#include "Node.h"

class Tree
{
  public:
    Tree();
    ~Tree();
    void init();
    void readFile(const string fn);
    void pruneSingles(Node *n); //recursive
    void readSubfams(string fn);
    void setSubfamNodes(Node *n); //recursive
    string getName();
    void setName(const string n);
    Node* getRoot();
    void setRoot(Node *r);
    string toNewickString();
    string newickStringForNode(Node *n);
    string nextNodeName();
    stack<string> leafNames(Node *n);
    void leafNodeNames(Node *n, stack<string> *names);
    vector<Node*> leaves(Node *n);
    void nodeLeaves(Node *n, vector<Node*> *leaves);
    void setLevelFrom(Node *n);
    int maxLevel();
    int getMaxNodeLevel(Node *n, int maxLevel);

  private:
    string name;
    Node *root;
    int nodeCount;
};

#endif

