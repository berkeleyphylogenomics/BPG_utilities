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

#ifndef TREE_SVG_H
#define TREE_SVG_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#include "SVG.h"
#include "Tree.h"

#define LABEL_WIDTH 50
#define NODE_COLOR  "gray"
#define SUBFAM_NODE_COLOR "red"

class TreeSVG
{
  public:
    TreeSVG();
    TreeSVG(int svgWide);
    TreeSVG(int svgWide, bool webMode);
    void  init(int svgWide, bool webMode);
    string toString(Tree *tree);
    void setNodeSVG(SVG *svg, Node *n);
    void setCoordinatesFrom(Node *n, int scaleX, int scaleY, int width, int &height);
    void initColors();
    void printImageMap(Tree *tree,string imfn);

  private:
    string *colors;
    int svgWidth;
    bool inWebMode;
};

#endif

