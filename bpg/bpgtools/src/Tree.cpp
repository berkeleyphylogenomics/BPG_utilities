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

#include "Tree.h"

Tree::Tree()
{
    init();
}

Tree::~Tree()
{
    delete root;
}

void Tree::init()
{
    nodeCount = 99;
    name      = "";
    root      = new Node();
}

void Tree::readFile(const string fn)
{
//    cout << "fn: " << fn << endl;
    //setName(fn);

    ifstream fin;
    fin.open(fn.c_str());

    if (! fin)
    {
        cerr << "ERROR: Could not open file: " << fn << endl;
        exit(1);
    }
    
    setName(fn);
    int  nodeCount = -1;
    int leafCount = -1;
    bool inBootstrap = false;
    string nodeName  = "";
    string bootstrap = "";
    Node *parent     = root;
    Node *child      = NULL;

    root->setIndex(++nodeCount);

    char ch;
    char ch2;

    while (fin.get(ch))
    {
        if (ch == '\n' || ch == ' ')
        {
            continue;
        }
        else if (ch == '(')
        {
            inBootstrap = false;

            child = new Node();
            child->setParent(parent);
            //child->setIndex(++nodeCount);
            parent->addChild(child);
            parent = child;
        }
        else if (ch == ')' || ch == ',')
        {
            inBootstrap = false;

            if (nodeName != "")
            {
                
//            parent->setBootstrap(bootstrap);
//            child = new Node();
//            child->setIndex(++leafCount);
                if(nodeName.size() > 3 && nodeName.substr(nodeName.size()-3,2) == "SF") //1 digit SF
                {
		//we have a labeled sequence
                    parent->setOtherSubfam(atoi(nodeName.substr(nodeName.size()-1).c_str()));
		//chop off family num
                    nodeName = nodeName.substr(0,nodeName.size()-4);
                }
                else if(nodeName.size() > 4 && nodeName.substr(nodeName.size()-4,2) == "SF") //2 digit SF
                {
		//we have a labeled sequence
                    parent->setOtherSubfam(atoi(nodeName.substr(nodeName.size()-2).c_str()));
		//chop off family num
                    nodeName = nodeName.substr(0,nodeName.size()-5);
                }
                else if(nodeName.size() > 5 && nodeName.substr(nodeName.size()-5,2) == "SF") //3 digit SF
                {
		//we have a labeled sequence
                    parent->setOtherSubfam(atoi(nodeName.substr(nodeName.size()-3).c_str()));
		//chop off family num
                    nodeName = nodeName.substr(0,nodeName.size()-6);
                }
                parent->setName(nodeName);
                parent->setIndex(++leafCount);
                parent->setLeaf(true);
            }
            else
            {
                parent->setIndex(++nodeCount);
            }
            //parent->setParent(parent);
            parent->setBootstrap(bootstrap);
            //parent->addChild(child);

            nodeName  = "";
            bootstrap = "";

            parent = parent->getParent();
            if (ch == ',')
            {
                child = new Node();
             //   child->setIndex(++leafCount);
                child->setParent(parent);
                parent->addChild(child);
                parent = child;
            }
        }
        else if (ch == ':')
        {
            inBootstrap = true;
        }
        else if (! inBootstrap)
        {
            // DB: WTF?!?  don't we want the ids to *not* change??
            ch2 = fin.peek();
            if((ch2 == ':' || ch2 == '\n' || ch2 == '\r') && ch == '_') //don't add extra _'s 
            {
                continue;
            }
            nodeName += ch;
        }
        else if (inBootstrap)
        {
            bootstrap += ch;
        }
    }

    fin.close();

    setLevelFrom(root);
}

/*Code added 7/12
 *pruneSingles() removes the extra node on singleton subfamilies
 */
void Tree::pruneSingles(Node *n) 
{
    if(n->isExternal())
    {
        return; //don't goof with external nodes
    }
    else if(n->getChildCount() != 1)
    {
        for (int i = 0; i < n->getChildCount(); i++)
        {
            pruneSingles(n->getChild(i));
        }
//        pruneSingles(n->getLeft());
//        pruneSingles(n->getRight());
    }
    else //must be singleton subfamily
    {
        n->getChild(0)->setParent(n->getParent());
        
        int index = n->getParent()->getChildIndex(n);   // get which slot we are are in in the parent.
        
        n->getParent()->setChild(n->getChild(0), index);

        n->removeChildren();
        delete n;
    }
}

/*marks certain nodes as the subfamily node
 *useful for drawing tree with colored dots for 
 *subfamily roots.
 *NOTE: Only useful after running readSubfams() or during
 *a kerf run!
 *NOTE: pruneSingles will alter the definition of a subfamily
 *before pruning, single member subfamilies have an extra node
 *which is the subfamily nod.
 *post pruning, that node disappears, so single subfamilies have
 *no such "subfamily node"
 *
 *pre-pruning is consistant with old versions of the code and makes
 *it consistant with gtree, etc
 *
 *pruning was added 7/05 to adapt the code for the web
 */
void Tree::setSubfamNodes(Node *n)
{
    //cout << "This node has " << n->getChildCount() << " children" << endl;
    
    for (int i = 0; i < n->getChildCount(); i++)
    {
        setSubfamNodes(n->getChild(i));
    }
    bool isSubfamNodeofSingle = ( n->getChildCount() < 2 ) && ( n->isInternal() );
    bool bothChildrenSameSubfam = ( n->getChildCount() > 1 );
    if (bothChildrenSameSubfam)
    {
        int sf = n->getChild(0)->getSubfam();
        for (int i = 0; i < n->getChildCount(); i++)
        {
            if ( i > 0 && n->getChild(i)->getSubfam() != sf)
            {
                bothChildrenSameSubfam = false;
                continue;
            }
        }
    }
    
//    ) && ( n->getLeft()->getSubfam() == n->getRight()->getSubfam() );
    bool parentDiffSubfam = (n->getParent() == NULL) || (  n->getParent()->getSubfam() != n->getSubfam() );
    //Singletons are single-member subfamilies with no subfamily node (eg post-pruned single member families)
    bool isSingleton = n->isExternal() && (n->getParent() != NULL && n->getParent()->getSubfam() != n->getSubfam());
    if(isSubfamNodeofSingle || (bothChildrenSameSubfam && parentDiffSubfam) )
    {
        n->setSubfamNode(true);
    }
    else
    {
        n->setSubfamNode(false);
    }
    if(isSingleton) 
    {
        n->setSingleton(true);
    } 
    else 
    {
        n->setSingleton(false);
    }
}

/*Method below added by Ryan Ritterson 7/5/05
 *readSubfamilies reads subfamilies from a .subfam file
 *and sets the subfamily of the coressponding node
 *note that the method in node will also set children
 *to same subfamily */
void Tree::readSubfams(string fn)
{
    ifstream fin;
    fin.open(fn.c_str());

    if (! fin)
    {
        cerr << "ERROR: Could not open subfam file: " << fn << endl;
        exit(1);
    }

    char ch;
    char prevch;
    int subfamilyNum = 1;
    string sequence = "";
    bool onSubfamily = false;
    bool onName = false;
    bool isFirstSeq = false;
    bool isFirstFam = true;
    Node* ancestor = NULL;

    while(fin.get(ch)) 
    {
        if((ch == '\n' || ch == ' ') && onName == false)
        {
            prevch = ch;
            continue;
        }
        else if(ch == '\r' || ( (ch == '(' || ch == ')') && onName == true) )
        {
            prevch = ch;
	    continue;
	}
        else if(ch == '%')
        {
	  //new subfamily
            if(!isFirstFam) 
	    {
                ancestor->setSubfam(subfamilyNum);
	    }
            else
            {
                isFirstFam = false;
            }
            onSubfamily = true;
            subfamilyNum++;
            isFirstSeq = true;
        }
        else if(ch == '>' && prevch == '\n')
        {
          //new sequence
            onName = true;
            onSubfamily = false;
            sequence = "";
        }
        else if(ch == ':' || ch == ',' || (ch == '\n' && onName == true)) 
        {
          //fell off sequence name
            onName = false;
          //now, search tree and find index of node with name
            Node *n = getRoot();
            Node *n1 = NULL;
            n->findNodeWithName(n1, sequence);
            if (n1 != NULL)
                n = n1;
            else
                cout << "Could not find node with name " << sequence << endl;
            
            if(isFirstSeq)
            {
                ancestor = n;
                isFirstSeq = false;
            }
            else 
            {
                ancestor = ancestor->getCommonAncestor(n);
            }
        }
        else if(onName) 
        {
            if(ch == '/' || ch == ' ')
            {
                ch = '_'; //bete does this, need it to match
            }
            sequence += ch;
        }
        prevch = ch;
    }
    ancestor->setSubfam(subfamilyNum);
}
	      
string Tree::getName()
{
    return name;
}

void Tree::setName(const string n)
{
//  cout << "n: " << n << endl;
    name = n;
}

Node* Tree::getRoot()
{
    return root;
}

void Tree::setRoot(Node *r)
{
    root = r;
}

string Tree::toNewickString()
{
    //return "(" + newickStringForNode(root) + ");";
    return newickStringForNode(root) + ";";
}

string Tree::newickStringForNode(Node *n)
{
    if (n->isExternal())
    {
        return (n->getBootstrap() == "") ? n->getName() : n->getName() + ":" + n->getBootstrap();
    }

    ostringstream ret;
    ret << "(";
    for (int i = 0; i < n->getChildCount(); i++)
    {
        if (i > 0)
            ret << ",";

        ret << newickStringForNode(n->getChild(i));
    }
    ret << ")";

    // DB: I don't have that pct_id crap in there yet (if ever).
    if (n->getBootstrap() != "" && n->getParent() != NULL)
        ret << ":" << n->getBootstrap();

    return ret.str();
}

string Tree::nextNodeName()
{
    nodeCount++;
    return "N";
}

stack<string> Tree::leafNames(Node *n)
{
    // Traverse the tree in inorder

    stack<string> names;

    leafNodeNames(n, &names);

    return names;
}

void Tree::leafNodeNames(Node *n, stack<string> *names)
{
//    stack<Node*> s;

//    for ( ; ; )
//    {
    if (n->isExternal())
        names->push(n->getName());
    else
    {
        for (int i = 0; i < n->getChildCount(); i++)
        {
            leafNodeNames(n->getChild(i), names);
        }
    }
}

vector<Node*> Tree::leaves(Node *n)
{
    vector<Node*> leaves;
    
    nodeLeaves(n, &leaves);

    return leaves;
}

void Tree::nodeLeaves(Node *n, vector<Node*> *leaves)
{
    if (n->isExternal())
    {
        leaves->push_back(n);
    }
    else
    {
        for (int i = 0; i < n->getChildCount(); i++)
        {
            nodeLeaves(n->getChild(i), leaves);
        }
    }
}

void Tree::setLevelFrom(Node *n)
{
    n->setLevel(0);
}

int Tree::maxLevel()
{
    // Traverse the tree in inorder

    int maxLevel = 0;

    Node *n = getRoot();
    return getMaxNodeLevel(n, maxLevel);
}

int Tree::getMaxNodeLevel(Node *n, int maxLevel)
{
    int nodeLevel = n->getLevel();
    
    maxLevel = nodeLevel > maxLevel ? nodeLevel : maxLevel;
    
    for (int i = 0; i < n->getChildCount(); i++)
    {
        maxLevel = getMaxNodeLevel(n->getChild(i), maxLevel);
    }
    
    return maxLevel;
}

