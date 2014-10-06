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

#include <math.h>

#include "TreeSVG.h"

TreeSVG::TreeSVG()
{
  init(320,false);
}

TreeSVG::TreeSVG(int svgWide)
{
  init(svgWide,false);
}

TreeSVG::TreeSVG(int svgWide, bool webMode)
{
  init(svgWide, webMode);
}

void TreeSVG::init(int svgWide, bool webMode)
{
  initColors();
  svgWidth = svgWide;
  inWebMode = webMode;
}

string TreeSVG::toString(Tree *tree)
{
    int fontSize = 8;
    int scaleY   = 12;
    //int svgWidth = 320;
    int svgHeight = tree->leafNames(tree->getRoot()).size() * scaleY + 2 * scaleY;
    int scaleX = (svgWidth - LABEL_WIDTH) / (tree->maxLevel() + 2);

    int height = 0;
    
    setCoordinatesFrom(tree->getRoot(), scaleX, scaleY, svgWidth, height);

    SVG svg;

    svg.setWidth(svgWidth);
    svg.setHeight(svgHeight);
    svg.setScaleX(scaleX);
    svg.setScaleY(scaleY);

    string stroke;
    string fontColor;

    stack<string> names;

    Node *n = tree->getRoot();
    stack<Node*> s;
    tree->setSubfamNodes(tree->getRoot()); //ensure subfamnodes are marked for proper coloring

    //    int childIndex = 0;
    for ( ; ; )
    {
        if (n != NULL)
        {
            //childIndex = 0;
            s.push(n);
            n = n->getChild(0);
        }
        else //n is null now, so pop off of stack
        {
            if (s.size() == 0)
            {
                break;
            }
          
            n = s.top();
            s.pop();
                
            // Begin visit
          
            if (n->isExternal())
            {
                float labelY = n->getY() + fontSize / 2;
                float labelX = n->getX() + 4;
              
                svg.addComment("Leaf: " + n->getName());
                //if we want to color text, find out what color to use
                if (n->getOtherSubfam() >= 0)
                {
                    int theColor = n->getOtherSubfam();
                
                    while (theColor > 31)
                    {
                        theColor -= 31;
                    }
                
                    fontColor = colors[theColor];
                }
                else
                {
                    fontColor = "black";
                }
                svg.addText(labelX, labelY, n->getName(),fontColor);
            }

            if (n->getParent() != NULL)
            {
              
                svg.addComment("Line left to parent x");
                if (n->getSubfam() >= 0)
                {
             
                    //Node *sib = n->nextSibling();
                    
                    if (!n->isSubfamNode() || (  n->isExternal() ) && 
                           ( n->nextSibling() != NULL && n->nextSibling()->getSubfam() != n->getSubfam() ) )
                    {
                  
                        int theColor = n->getSubfam();
                  
                        while (theColor > 31)
                        {
                            theColor -= 31;
                        }
                        stroke = colors[theColor];
                    }
                    else 
                    {
                        stroke = "black";
                    }

                }
                else
                {
                    stroke = "black";
                }
              
                svg.addPath(n->getX(), n->getY(), n->getParent()->getX(), n->getY(), stroke);
              
                svg.addComment("Line Up/down to parent");
              
                svg.addPath(n->getParent()->getX(), n->getY(),
                            n->getParent()->getX(), n->getParent()->getY(), stroke);
            }
            else
            {
                svg.addComment("Root");
              
                stroke = "black";
              
                svg.addPath(n->getX(), n->getY(), (n->getX() - scaleX),
                            n->getY(), stroke);
            }

            // End visit

            n = n->nextSibling();
            /*
            Node *parent = n->getParent();
            if (parent != NULL)
            {
                int childIndex = parent->getChildIndex(n);
                n = n->getChild(childIndex + 1);
            }
 //          n = n->getChild(++childIndex);
            */
        }
    }
    setNodeSVG(&svg, tree->getRoot());
  
    return svg.toString();
}

/**
 **prints an HTML image map suitable for direct insertion into an HTML file
 **for the SVG
 **/

void TreeSVG::printImageMap(Tree *tree,string imfn) //image map file name
{
    int fontSize = 8;
    int scaleY   = 12;
    //int svgWidth = 320;
    int scaleX = (svgWidth - LABEL_WIDTH) / (tree->maxLevel() + 2);

    ofstream imagemap;
    imagemap.open(imfn.c_str());
    float subFamPos[1000][4]; //int[subfamNum][xmin,ymin,xmax,ymax]
    Node* subFamNodes[1000];

    for(int a=0; a<= 999; a++) 
    {
        subFamPos[a][0] = 100000.0;
        subFamPos[a][1] = 100000.0;
        subFamPos[a][2] = 0.0;
        subFamPos[a][3] = 0.0;
        subFamNodes[a] = NULL;
    }

    Node* n = tree->getRoot();
    stack<Node*> s;
//    int index = 0;
    //traverse tree to get subfamilies and coordinates
    while(true)
    {
        if (n != NULL)
        {
            //index = 0;
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
            int subFam = n->getSubfam();

            if (subFam >= 0)
            {
                if(n->isSubfamNode() || n->isSingleton())
                {
                    subFamNodes[subFam] = n;
                }

                float x = n->getX();
                float y = n->getY();

                if (x < subFamPos[subFam][0])
                {
                    subFamPos[subFam][0] = x;
                }

                if (x > subFamPos[subFam][2])
                {
                    subFamPos[subFam][2] = x;
                }

                if (y < subFamPos[subFam][1])
                {
                    subFamPos[subFam][1] = y;
                }

                if (y > subFamPos[subFam][3])
                {
                    subFamPos[subFam][3] = y;
                }
            }
      
            // End visit
            n = n->getSibling();
            /*
            Node *parent = n->getParent();
            if (parent != NULL)
            {
                int childIndex = parent->getChildIndex(n);
                n = n->getChild(childIndex + 1);
            } */         
            // n = n->getChild(++index);
        }
    }

    //now, loop over array and print all image map values.
   for(int a = 0; a<=999; a++)
   {
        if (subFamPos[a][0] == 100000)
        {
            continue; //no values
        }
        else
        {
            if (subFamNodes[a]->isSingleton()) //singleton subfamily
            {
                subFamPos[a][0] = subFamNodes[a]->getParent()->getX();
                subFamPos[a][1] -= 2;
                subFamPos[a][3] += 2; //slightly larger than node circles, center on line
            }

            imagemap << a << ",coords=\"" << (subFamPos[a][0]) << ",";
            imagemap << subFamPos[a][1] << "," << subFamPos[a][2];
            imagemap << "," << subFamPos[a][3] << "\"\n";
        }
     }

    imagemap.close();
}
  
void TreeSVG::setNodeSVG(SVG *svg, Node *n)
{
    if (n == NULL)
    {
        return;
    }

    if (n->isInternal())
    {
        string color;

        if (n->isSubfamNode())
        {
            color = SUBFAM_NODE_COLOR;
        }
        else
        {
            color = NODE_COLOR;
        }

        svg->addCircle(n->getX(), n->getY(), 3, color, color);

        if (! inWebMode) 
        {
            svg->addText(n->getX() + 6, n->getY() + 4, n->getBootstrap(), 6);
        }
    }

    for (int i = 0; i < n->getChildCount(); i++)
    {
        setNodeSVG(svg, n->getChild(i));
    }
}

void TreeSVG::setCoordinatesFrom(Node *n, int scaleX, int scaleY, int width, int &height)
{
    if (n == NULL)
    {
        return;
    }

    if (n->isInternal())
    {
        for (int i = 0; i < n->getChildCount(); i++)
        {
            setCoordinatesFrom(n->getChild(i), scaleX, scaleY, width, height);
        }
        
        float x;
        float y;

        if (n->getChildCount() == 1)
        {
            x = n->getChild(0)->getX();
            y = n->getChild(0)->getY();
        }
        else
        {
            x = 1e10;
            float minY = 1e10;
            float maxY = -1;
            
            for (int i = 0; i < n->getChildCount(); i++)
            {
                x = fmin(x, n->getChild(i)->getX());
                minY = fmin(minY, n->getChild(i)->getY());
                maxY = fmax(maxY, n->getChild(i)->getY());
            }
            y = (minY + maxY) /2;
        }
        /*
        if (n->hasRight() && n->hasLeft())
        {
            x = fmin(n->getRight()->getX(), n->getLeft()->getX());
            y = (n->getRight()->getY() + n->getLeft()->getY()) / 2;
        }
        else if (n->hasRight())
        {
            x = n->getRight()->getX();
            y = n->getRight()->getY();
        }
        else if (n->hasLeft())
        {
            x = n->getLeft()->getX();
            y = n->getLeft()->getY();
        }

        */
        x -= scaleX;
        n->setX(x);
        n->setY(y);
    }
    else
    {
        n->setX(width - LABEL_WIDTH - scaleX);
        height += scaleY;
        n->setY(height);
    }
}

void TreeSVG::initColors()
{
    colors = new string[32];

    colors[0] = "crimson";
    colors[1] = "forestgreen";
    colors[2] = "cornflowerblue";
    colors[3] = "yellow";
    colors[4] = "darkorange";
    colors[5] = "blueviolet";
    colors[6] = "chartreuse";
    colors[7] = "darkkhaki";
    colors[8] = "gold";
    colors[9] = "burlywood";
    colors[10] = "coral";
    colors[11] = "deeppink";
    colors[12] = "dodgerblue";
    colors[13] = "cyan";
    colors[14] = "darkorange";
    colors[15] = "indigo";
    colors[16] = "lightskyblue";
    colors[17] = "mediumpurple";
    colors[18] = "mediumspringgreen";
    colors[19] = "olivedrab";
    colors[20] = "orchid";
    colors[21] = "salmon";
    colors[22] = "slateblue";
    colors[23] = "tomato";
    colors[24] = "yellowgreen";
    colors[25] = "darkcyan";
    colors[26] = "sienna";
    colors[27] = "thistle";
    colors[28] = "peru";
    colors[29] = "palegreen";
    colors[30] = "orange";
    colors[31] = "saddlebrown";
}

