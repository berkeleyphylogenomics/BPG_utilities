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

#include "Sequence.h"

Sequence::Sequence()
{
    init();
}

void Sequence::init()
{
    id = "";
    definition = "";
    sequence = "";
}

string Sequence::getDefinition()
{
    return definition;
}

void Sequence::setDefinition(string newDefinition)
{
    definition = newDefinition;
    setIdFromDefinition();
}

string Sequence::getId()
{
    return id;
}

void Sequence::setId(string newId)
{
    id = newId;
}

string Sequence::getSequence()
{
    return sequence;
}

void Sequence::setSequence(string newSequence)
{
    sequence = newSequence;
}

int Sequence::length()
{
    return sequence.length();
}

string Sequence::getResidue(int index)
{
    return sequence.substr(index, 1);
}

string Sequence::toFasta()
{
    return definition + "\n" + sequence;
}

void Sequence::setIdFromDefinition()
{
    string tmpid = definition;

    if ( tmpid.substr(0,1) == ">")
    {
        tmpid = tmpid.substr(1, tmpid.length());
    }
    
    string::size_type pos = tmpid.find(" ", 0);

    if (pos != string::npos)
    {
        tmpid = tmpid.substr(0, pos);
    }
    
// DB: WTF?!?!??
    // It seems like the code *intentionally* removes characters from the ids!  Dumb!
    /*
    pos = definition.find("|", 0);

    if (pos != string::npos)
    {
        tmpid = tmpid.substr(pos, tmpid.length());
    }
    */

    id = tmpid;
}

