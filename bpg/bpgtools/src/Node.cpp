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

#include "Node.h"

Node::Node()
{
    init();
}

Node::~Node()
{
    /*if (isInternal())
    {
        cout << "Deleting internal node " << index << endl;
    }
    else
    {
        cout << "Deleting leaf " << index << endl;
    }*/
    if (children != NULL)
    {
        for (int i = 0; i< childCount; i++)
        {
            Node *n = getChild(i);
            delete n;
        }
        delete [] children;
    }
}

void Node::init()
{
    name      = "";
//    right     = NULL;
//    left      = NULL;
    children  = NULL;
    parent    = NULL;
    childCount = 0;
    bootstrap = "100";
    pid       = 100.0;
    leaf      = false;
    index     = -1;
    level     = 0;
    subfam    = -1;
    otherSubfam = -1;
    subfamNode = false;
    singleton = false;
}

bool Node::isBinary()
{
    return (childCount == 2);
}

/*bool Node::hasRight()
{
    return right != NULL;
}

bool Node::hasLeft()
{
    return left != NULL;
}
*/
bool Node::isInternal()
{
    //return hasRight() || hasLeft();
    return !leaf;
}

bool Node::isExternal()
{
    //return ! (hasRight() || hasLeft());
    return leaf;
}

void Node::setLeaf(bool l)
{
    leaf = l;
}

string Node::getName()
{
    return name;
}

void Node::setName(string n)
{
    name = n;
}

Node* Node::getParent()
{
    return parent;
}

void Node::setParent(Node *n)
{
    parent = n;
}
/*
Node* Node::getRight()
{
    return right;
}

void Node::setRight(Node *n)
{
    right = n;
}

Node* Node::getLeft()
{
    return left;
}

void Node::setLeft(Node *n)
{
    left = n;
}
*/
string Node::getBootstrap()
{
    return bootstrap;
}

int Node::getBootstrapInt()
{
    return atoi(bootstrap.c_str());
}

void Node::setBootstrap(string b)
{
    int dot = b.find_first_of(".", 0);

    if (dot != string::npos)
    {
        bootstrap = b.substr(0, dot).c_str();
    }
    else
    {
        bootstrap = b.c_str();
    }
}

void Node::removeChildren()
{
//    right = NULL;
//    left  = NULL;
/*    if (children != NULL)
    {
        for (int i = 0; i < childCount; i++)
        {
            children[i]->removeChildren();
        }
        delete[] children;
    }
    */
    delete [] children;
    children = NULL;
    childCount = 0;
}
/*
void Node::removeLeft()
{
    left = NULL;
}

void Node::removeRight()
{
    right = NULL;
}
*/

bool Node::findNodeWithName(Node *n, string seqName)
{
    if (isExternal() && seqName == name)
    {
        n = this;
        return true;
    }
    else
    {
        for (int i = 0; i < childCount; i++)
        {
            if (getChild(i)->findNodeWithName(n, seqName))
            {
                return true;
            }
        }
    }
    return false;
}

void Node::addChild(Node *n)
{
    if (childCount == 0)
    {
        children = new Node*[1];
        children[0] = n;
    }
    else
    {
        // allocate another slot for new child
        Node **tmp = children;
        children = new Node*[childCount+1];
        memcpy(children, tmp, childCount * sizeof(Node*));
        children[childCount] = n;
        delete tmp;
    }
    childCount++;
    /*
    if (! hasLeft())
    {
        setLeft(n);
    }
    else if (! hasRight())
    {
        setRight(n);
    }
    else
    {
        cerr << "ERROR: Attempt to multifurcate tree!" << endl;
        exit(1);
    }*/
}

float Node::getPercentId()
{
    return pid;
}

void Node::setPercentId(float p)
{
    pid = p;
}

int Node::getIndex()
{
    return index;
}

void Node::setIndex(int i)
{
    index = i;
}

int Node::getChildCount()
{
    return childCount;
}

void Node::setChild(Node *n, int i)
{
    if (childCount > i && i >= 0)
        children[i] = n;
}

Node *Node::getChild(int i)
{
    if (i < childCount && i >= 0)
    {
        return children[i];
    }
    return NULL;
}

int Node::getChildIndex(Node *n)
{
    for (int i = 0; i < childCount; i++)
    {
        if (children[i] == n )
            return i;
    }
    return -1;
}

Node *Node::nextSibling()
{
    if (parent != NULL)
    {
        int myIndex = parent->getChildIndex(this);
        return parent->getChild(myIndex+1);
    }
    return NULL;
}

Node *Node::getSibling()
{
    if (parent != NULL)
    {
        for (int i = 0; i < parent->getChildCount(); i++)
        {
            if (parent->getChild(i) != this)
            {
                return parent->getChild(i);
            }
        }
        
    /*
        if (parent->getRight() == this)
        {
            return parent->getLeft();
        }
        else
        {
            return parent->getRight();
        }
        */
    }
    return NULL;
}

int Node::getLevel()
{
    return level;
}

void Node::setLevel(int l)
{
    level = l;

    for (int i = 0; i < childCount; i++)
    {
        getChild(i)->setLevel(l + 1);
    }
    /*
    if (hasRight())
    {
        getRight()->setLevel(l + 1);
    }

    if (hasLeft())
    {
        getLeft()->setLevel(l + 1);
    }*/
}

/*Method coded by Ryan Ritterson
 *returns nearest common ancestor node
 *of this node and input node */
Node* Node::getCommonAncestor(Node* b)
{
  Node* a = this;
  
  while(a != b)
    {
      if( a->getLevel() > b->getLevel() )
	{
	  //a farther down, so move up
	  a = a->getParent();
	}
      else if ( a -> getLevel() < b-> getLevel() ) 
	{
	  //b farther down, so move up
	  b = b->getParent();
	}
      else 
	{
	  //at same level, not equal, so move both up
	  a = a->getParent();
	  b = b->getParent();
	}
    }
  return b;
}

float Node::getX()
{
    return xCoordinate;
}

void Node::setX(float x)
{
    xCoordinate = x;
}

float Node::getY()
{
    return yCoordinate;
}

void Node::setY(float n)
{
    yCoordinate = n;
}

int Node::getSubfam()
{
    return subfam;
}

void Node::setSubfam(int s)
{
    subfam = s;

    for (int i = 0; i < childCount; i++)
    {
        getChild(i)->setSubfam(s);
    }
    /*
    if (hasLeft())
    {
        getLeft()->setSubfam(s);
    }

    if (hasRight())
    {
        getRight()->setSubfam(s);
    }*/
}

int Node::getOtherSubfam()
{
  return otherSubfam;
}

//note that this function doesn't tamper with children.
//is is only intended to be run on leaf nodes.
void Node::setOtherSubfam(int o) 
{
  otherSubfam = o;
}

bool Node::isSubfamNode()
{
  return subfamNode;
}

void Node::setSubfamNode(bool b)
{
  subfamNode = b;
}

bool Node::isSingleton()
{
  return singleton;
}

void Node::setSingleton(bool b)
{
  singleton = b;
}
