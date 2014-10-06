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

#ifndef SVG_H
#define SVG_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

class SVG
{
  public:
    SVG();
    void   init();
    string toString();
    void   saveToFile(string fn);
    string getHeader();
    string getFooter();
    int    getWidth();
    void   setWidth(int w);
    int    getHeight();
    void   setHeight(int h);
    int    getScaleX();
    void   setScaleX(int x);
    int    getScaleY();
    void   setScaleY(int y);
    string getFontFamily();
    void   setFontFamily(string family);
    void   addText(float x, float y, string text);
    void   addText(float x, float y, string text, int size);
    void   addText(float x, float y, string text, string color);
    void   addText(float x, float y, string text, string color, int size);
    void   addPath(float origx, float origy, float destx, float desty);
    void   addPath(float origx, float origy, float destx, float desty,
        string stroke);
    void   addPath(float origx, float origy, float destx, float desty,
        string stroke, int strokeWidth);
    void   addCircle(float x, float y, int radius, string fill, string stroke);
    void   addComment(string c);

  private:
    int width;
    int height;
    int scaleX;
    int scaleY;
    int fontSize;
    string fontFamily;
    string fontWeight;
    string fontStyle;
    string fontFill;
    ostringstream svg;
};

#endif

