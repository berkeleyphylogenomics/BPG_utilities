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

#include "SVG.h"

SVG::SVG()
{
    init();
}

void SVG::init()
{
    width    = 640;
    height   = 480;
    scaleX   = 1;
    scaleY   = 1;
    fontSize = 8;
    fontFamily = "Courier";
}

string SVG::toString()
{
    return getHeader() + svg.str() + getFooter();
}

void SVG::saveToFile(string fn)
{
    ofstream svgfile;
    svgfile.open(fn.c_str());
    svgfile << toString();
    svgfile.close();
}

string SVG::getHeader()
{
    ostringstream header;

    header << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl << endl;
    header << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20010904//EN\"" << endl;
    header << "\"http://www.w3.org/TR/2000/REC-SVG-20010904/DTD/";
    header << "svg10.dtd\">" << endl << endl;

    header << "<svg width=\"" << width << "\" height=\"" << height <<
        "\" viewBox=\"0 0 " << width << " " << height << "\">" << endl;
    header << "<title>Subfamily Tree</title>" << endl;
    header << "<desc>Subfamily Tree</desc>" << endl;

    header << "<g style=\"fill:black; stroke:black; stroke-width:1\">" << endl;

    return header.str();
}

string SVG::getFooter()
{
    return "</g>\n</svg>\n";
}


void SVG::addText(float x, float y, string text)
{
  string black = "#000000";
  addText(x, y, text, black, fontSize);
}

void SVG::addText(float x, float y, string text, int size) 
{
  string black = "#000000";
  addText(x,y,text,black,size);
}

void SVG::addText(float x, float y, string text, string color)
{
  addText(x,y,text,color,fontSize);
}

void SVG::addText(float x, float y, string text, string color, int size)
{
  svg << "<text x=\"" << x << "\" y=\"" << y << "\" " <<
    "style=\"font-family:" << fontFamily << "; font-weight:normal; " <<
    "font-style:normal; font-size:" << fontSize <<
    "pt; fill:" << color << "; " << "stroke:" << color <<
    "; stroke-width:0;\">" <<
    text << "</text>" << endl;
}

void SVG::addPath(float origx, float origy, float destx, float desty)
{
    addPath(origx, origy, destx, desty, "black", 2);
}

void SVG::addPath(float origx, float origy, float destx, float desty,
    string stroke)
{
    addPath(origx, origy, destx, desty, stroke, 2);
}

void SVG::addPath(float origx, float origy, float destx, float desty,
    string stroke, int strokeWidth)
{
    svg << "<path d=\"M" << origx << " " << origy << " L" << destx <<
        " " << desty << "\" stroke=\"" << stroke << "\" stroke-width=\"" <<
        strokeWidth << "\"/>" << endl;
}

void SVG::addCircle(float x, float y, int radius, string fill, string stroke)
{
    svg << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << radius <<
        "\" fill=\"" << fill << "\" stroke=\"" << stroke << "\"/>" << endl;
}

void SVG::addComment(string c)
{
    svg << "<!-- " << c << " -->" << endl;
}

int SVG::getWidth()
{
    return width;
}

void SVG::setWidth(int w)
{
    width = w;
}

int SVG::getHeight()
{
    return height;
}

void SVG::setHeight(int h)
{
    height = h;
}

int SVG::getScaleX()
{
    return scaleX;
}

void SVG::setScaleX(int x)
{
    scaleX = x;
}

int SVG::getScaleY()
{
    return scaleY;
}

void SVG::setScaleY(int y)
{
    scaleY = y;
}

string SVG::getFontFamily()
{
    return fontFamily;
}

void SVG::setFontFamily(string name)
{
    fontFamily = name;
}

