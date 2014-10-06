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

#ifndef ALIGNMENT_STATISTICS_H
#define ALIGNMENT_STATISTICS_H

#include <string>
#include <iostream>
using namespace std;

#include "Fasta.h"
#include "Blosum62.h"

class AlignmentStatistics
{
  public:
    AlignmentStatistics();
    void init();

    float pairPercentId(Sequence *s, Sequence *t);
    int compute(Fasta f);
    int checkLengths(Fasta f);
    int global(Fasta f);
    int core(Fasta f);

    int   numberOfSequences();
    int   alignmentLength();
    int   alignedColumns();

    string shortestSequenceId();
    string shortestAlignedRegionSequenceId();

    string longestSequenceId();
    string longestAlignedRegionSequenceId();

    float averageLength();
    float averageAlignedRegionLength();

    int   shortestLength();
    int   shortestAlignedRegionLength();

    int   longestLength();
    int   longestAlignedRegionLength();

    float meanPercentIdentity();

    float globalPercentIdentity();
    float approximateMeanPercentIdentity();

    float  minimumPercentIdentity();
    string minimumPercentIdentityIds();

    float  maximumPercentIdentity();
    string maximumPercentIdentityIds();
    
    float percentGaps();
    float percentColumnsNegativeBlosum62();
    int   columnsNegativeBlosum62();

    int   maximumLengthDifference();
    int   maximumAlignedRegionLengthDifference();

    float percentMaximumLengthDifference();
    float percentAlignedRegionMaximumLengthDifference();

    float meanAlignedToCore();
    int   coreLength();
    int   firstCoreResidue();
    int   lastCoreResidue();
    float meanCoreGapPercentage();

    void allowGaps(bool b);
    bool allowingGaps();

  private:
    float meanpi;
    float minpi;

    float globalpi;
    float approxmeanpi;

    float meanCore;
    float meanCoreGaps;
    int firstCoreRes;
    int lastCoreRes;
    int lengthCore;
    
    string minpiids;
    float maxpi;
    string maxpiids;
    float pgaps;
    float pnb62;
    int   cnb62;
    int   maxld;
    float pmaxld;
    Blosum62 blosum62;

    int  ncols;
    int  nseqs;

    int sumSlen;
    int sumSAlignedRegionLen;

    string shortest;
    string shortestAlignedRegion;

    string longest;
    string longestAlignedRegion;

    int shortestLen;
    int shortestAlignedRegionLen;

    int longestLen;
    int longestAlignedRegionLen;

    bool letGaps;
    int alignedRegionLength;
};

#endif

