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

#include <fstream>
#include <iostream>

#include "Fasta.h"

Fasta::Fasta()
{
    init();
}

Fasta::~Fasta()
{
}

void Fasta::init()
{
    sequenceCount = 0;
    sequences = new Sequence[MAX_SEQS];
}

int Fasta::readFile(string filename)
{
    ifstream ifs( filename.c_str() );
    if ( ! ifs )
    {
        cerr << "Could not open file " << filename << endl;
        return( 1 );
    }

    string def;
    string seq;
    string line;

    while (getline(ifs, line))
    {
        if (line.length() == 0)
        {
            continue;
        }

        if (line.substr(0, 1) == ">")
        {
            if (seq.length() > 0)
            {
                Sequence s;
                s.setDefinition(def);
                s.setSequence(seq); 

                addSequence(s);

                seq.erase(0, seq.length());
            }

            def.erase(0, def.length());
            def = line;
        }
        else
        {
            seq += line;
        }

        line.erase(0, line.length());
    }

    Sequence s;
    s.setDefinition(def);
    s.setSequence(seq);
    addSequence(s);
}

int Fasta::count()
{
    return sequenceCount;
}

Sequence Fasta::getSequence(int index)
{
    if (index < sequenceCount)
    {
        return sequences[index];
    }
}

Sequence Fasta::getSequence(string definitionOrId)
{
    string::size_type pos;

    for (int i = 0; i < sequenceCount; i++)
    {
        Sequence s = sequences[i];

        pos = s.getId().find(definitionOrId, 0);

        if (pos != string::npos)
        {
            return s;
        }
    }

    Sequence s;
    return s;
}

void Fasta::addSequence(Sequence sequence)
{
    if (sequenceCount < MAX_SEQS - 1)
    {
        sequenceCount++;

        sequences[sequenceCount - 1] = sequence;
    }
    else
    {
        cerr << "ERROR: Attempt to add too many sequences" << endl;
    }
}
/*
void Fasta::removeSequence(int index)
{
    if (index < sequenceCount && index < MAX_SEQS - 1)
    {
        seqCount--;

        for (int i = index; i < seqCount; i++)
        {
            sequences[i] = sequence[i + 1];
        }
    }
    else
    {
        cerr << "ERROR: Attempt to remove bad index: " << index << endl;
    }
}
*/
string Fasta::toFasta()
{
    string fasta;

    for (int i = 0; i < sequenceCount; i++)
    {
        Sequence s = getSequence(i);
        fasta += s.toFasta() + "\n";
    }

    return fasta;
}

