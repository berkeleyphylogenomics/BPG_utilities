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

#ifndef NODE_H
#define NODE_H

#include <string>
#include <iostream>
using namespace std;

class Node
{
  public:
    Node();
    ~Node();
    void init();
    bool isBinary();
/*    bool hasRight();
    bool hasLeft(); */
    bool isInternal();
    bool isExternal();
    string getName();
    void setName(string n);
    Node *getParent();
    void setParent(Node *n);
/*    Node *getRight();
    void setRight(Node *n);
    Node *getLeft();
    void setLeft(Node *n); */
    string getBootstrap();
    int getBootstrapInt();
    void setBootstrap(string b);
    void addChild(Node *n);
    float getPercentId();
    void setPercentId(float p);
    void removeChildren();
//    void removeLeft();
//    void removeRight();
    int getIndex();
    void setIndex(int i);
    Node *getSibling();
    void setLevel(int l);
    int getLevel();
    Node* getCommonAncestor(Node* b);
    float getX();
    void setX(float x);
    float getY();
    void setY(float y);
    int getSubfam();
    void setSubfam(int s);
    void setOtherSubfam(int o);
    int getOtherSubfam();
    void setSubfamNode(bool b);
    bool isSubfamNode();
    void setSingleton(bool b);
    bool isSingleton();
// Duncan:
//    void removeChild(int childIndex);
//    Node *getSibling(int sibIndex);
    Node *getChild(int i);
    void setChild(Node* n, int childIndex);
    void setLeaf(bool l);
    int getChildCount();
    int getChildIndex(Node *n);
    bool findNodeWithName(Node *n, string seqName);
    Node *nextSibling();
    
  private:
    string name;
    string bootstrap;
//    Node   *left;
//    Node   *right;
    Node   **children;
    Node   *parent;
    int    childCount;
    float  pid;
    int    index;
    int    level;
    float  xCoordinate;
    float  yCoordinate;
    int    subfam;
    int    otherSubfam;
    bool   subfamNode;
    bool   singleton;
    bool   leaf;
};

#endif

