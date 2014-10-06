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

#include "SVG.h"
#include "AlignmentSVG.h"

AlignmentSVG::AlignmentSVG()
{
    init();
}

void AlignmentSVG::init()
{
}

string AlignmentSVG::toString(Fasta *fa)
{
    int rowCount = fa->count();
    int colCount = fa->getSequence(0).length();

    int fontSize = 8;
    int scaleX   = 8;
    int scaleY   = 12;
    int svgWidth = colCount * scaleX + 5;
    int svgHeight = rowCount * scaleY + 5;

    SVG svg;

    svg.setWidth(svgWidth);
    svg.setHeight(svgHeight);
    svg.setScaleX(scaleX);
    svg.setScaleY(scaleY);

    Sequence s;

    for (int row = 0; row < rowCount; row++)
    {
        s = fa->getSequence(row);

        svg.addText(5, row * scaleY + (fontSize * 2) , s.getSequence());
    }

    return svg.toString();
}

