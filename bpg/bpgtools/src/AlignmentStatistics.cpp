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

#include <string>
#include <iostream>
using namespace std;

#include "AlignmentStatistics.h"

AlignmentStatistics::AlignmentStatistics()
{
    init();
}

void AlignmentStatistics::init()
{
    minpi  = 0.0;

    globalpi = 0.0;
    approxmeanpi= 0.0;
    meanpi = 0.0;

    maxpi  = 0.0;
    pgaps  = 0.0;
    pnb62  = 0.0;
    cnb62  = 0;
    maxld  = 0;
    pmaxld = 1.0;
    letGaps = false; //doesn't count any column with gaps as aligned by default.

    meanCore = 0.0;
    firstCoreRes = 0;
    lastCoreRes = 0;
    lengthCore = 0;
    
}

float AlignmentStatistics::pairPercentId(Sequence *s, Sequence *t)
{
    if (letGaps)
    {
      // Multiple lowercase characters in one sequence may correspond to a gap 
      // in the other sequence

      int i = 0;
      int j = 0;
      float id = 0;
      float aligned = 0;
      while (i < s->length() || j < t->length() )
      {
        if (i < s->length()) 
        {
          char aas = s->getResidue(i).c_str()[0];
          if (islower(aas) || aas == '.')
          {
            // Skip past lowercase letters and gaps in s
            i++;
          }
          else
          {
            if (j < t->length()) {
              char aat = t->getResidue(j).c_str()[0];
              if (islower(aat) || aat == '.') {
                // Skip past lowercase letters and gaps in t
                j++;
              } else {
                if (isupper(aas) || isupper(aat)) {
                  aligned++;
                  if (aas == aat) {
                    id++;
                  }
                }
                i++;
                j++;
              }
            } else {
              // Have exhausted all characters in t
              if (isupper(aas)) {
                aligned++;
              }
              i++;
            }
          }
        } 
        else 
        {
          // Have exhausted all characters in s
          // Since we're in the body of the loop, j < t->length()
          char aat = t->getResidue(j).c_str()[0];
          // Skip past lowercase letters and gaps in t
          if (!islower(aat) && aat != '.') {
            if (isupper(aat)) 
            {
              aligned++;
            }
          }
          j++;
        }
      }
      return id / aligned * 100;
    }
    else
    {
        if (s->length() != t->length())
        {
            cerr << "Sequences must be of the same length!" <<endl;
            return -1;
        }

        int len = s->length();

        float id = 0;
        float aligned = 0;
        for (int i = 0; i < len; i++)
        {
            char aas = s->getResidue(i).c_str()[0];
            if (islower(aas) || aas == '.')
            {
                continue;
            }
            
            char aat = t->getResidue(i).c_str()[0];
            
            bool isAligned;
            if (!letGaps)
            {
                isAligned = isupper(aas) || isupper(aat);
            }
            else
            {
                isAligned = isupper(aas) && isupper(aat);
            }

            if (isAligned)
            {
                aligned++;
                if (aas == aat)
                    id++;
            }
        }
        return id / aligned * 100;
    }
    // Should never get here
    return -1.0;
}
             
int AlignmentStatistics::compute(Fasta f)
{
    ncols = checkLengths(f);

    if (ncols == 0)
    {
        cout << "Sequences must all be of the same length!" << endl;
        return 0;
    }

    Sequence s;
    Sequence t;

    nseqs = f.count();

    meanpi = 0.0;
    minpi  = 100.0;
    maxpi  = 0.0;
    pgaps  = 0.0;
    cnb62  = 0;
    pnb62 = 0.0;

    int colScores[ncols];
    //initialize column scores, variables
    for (int i = 0; i < ncols; i++)
    {
        colScores[i] = 0;
    }

    shortestLen = 99999;
    shortestAlignedRegionLen = 99999;
    longestLen = 0;
    longestAlignedRegionLen = 0;

    sumSlen = 0;
    sumSAlignedRegionLen = 0;

    int lowScoreCount = 0;
    int dashCount = 0;
    int dotCount = 0;
    int upperCount = 0;

    int sequencesCompared = 0;
    alignedRegionLength  = 0;

    // Loop over sequences.

    for (int i = 0; i < nseqs; i++)
    {
        s = f.getSequence(i);

        int sDashCount  = 0;
        int sDotCount   = 0;
        int sUpperCount = 0;

        char aa;
        
        // Count # dots, upper, dashes in each sequence.

        for (int j = 0; j < s.length(); j++)
        {
            aa = s.getResidue(j).c_str()[0];

            if (isupper(aa))
            {
                upperCount++;
                sUpperCount++;
            }
            else if (aa == '-')
            {
                dashCount++;
                sDashCount++;
            }
            else if (aa == '.')
            {
                dotCount++;
                sDotCount++;
            }
        }

        // Aligned region length (ignore inserts: dots and lowercase) -- do for
        // first sequence only.

        if (i == 0)
        {
            alignedRegionLength = upperCount + dashCount;
        }

        // Length - including inserts.  Sum for average.

        int slen = s.length() - sDashCount - sDotCount;
        sumSlen += slen;
        
        // Length - excluding inserts (i.e, aligned region only).  Sum for
        // average.

        int sAlignedRegionLen = sUpperCount;
        sumSAlignedRegionLen += sAlignedRegionLen;
        
        // Record this sequence as longest or shortest, if is such.

        if (slen < shortestLen)
        {
            shortestLen = slen;
            shortest = s.getId();
        }

        if (slen > longestLen)
        {
            longestLen = slen;
            longest = s.getId();
        }

        // Again for aligned columns.

        if (sAlignedRegionLen < shortestAlignedRegionLen)
        {
            shortestAlignedRegionLen = sAlignedRegionLen;
            shortestAlignedRegion = s.getId();
        }

        if (slen > longestAlignedRegionLen)
        {
            longestAlignedRegionLen = sAlignedRegionLen;
            longestAlignedRegion = s.getId();
        }

        // Get new sequence and compare to current -- 
        // loop only over those higher than current.

        for (int j = i + 1; j < nseqs; j++)
        {
            t = f.getSequence(j); 

            int matches = 0;
            int aligns  = 0;

            char schar;
            char tchar;

            for (int k = 0; k < ncols; k++)
            {
                schar = s.getResidue(k).c_str()[0];
                tchar = t.getResidue(k).c_str()[0];
                
                bool isAligned;
                
                // Ignore all gaps-- both need to be upper case.

                if(letGaps == true) 
                {
                    isAligned = isupper(schar) && isupper(tchar);
                }

                // One uppercase counts, even if the other is a gap, but they 
                // won't match.

                else if(letGaps == false) 
                {
                    isAligned = isupper(schar) || isupper(tchar);  
                }
                if (isAligned)
                {
                    aligns++;
                    
                    // If same character exactly, they match.

                    if (schar == tchar)
                    {
                        matches++;
                    }

                    int scor = blosum62.score(schar, tchar);

                    colScores[k] += blosum62.score(schar, tchar);
                }
            }

            
            sequencesCompared++;

            float thispi = (float) matches / aligns;
            meanpi += thispi;

            if (thispi < minpi)
            {
                minpi = thispi;
                minpiids = s.getId() + ", " + t.getId();
            }

            if (thispi > maxpi)
            {
                maxpi = thispi;
                maxpiids = s.getId() + ", " + t.getId();
            }
        }
    }

    /* If only two sequences, max = min. */

    if (maxpi == 0.0)
    {
        maxpi = minpi;
    }

    meanpi /= sequencesCompared;
    meanpi *= 100;

    minpi *= 100;
    maxpi *= 100;

    pgaps = (float) dashCount / (dashCount + upperCount);
    pgaps *= 100;

    for (int i = 0; i < ncols; i++)
    {
        float avg = (float) colScores[i] / nseqs;

        if (avg < 0)
        {
            cnb62++;
        }
    }

    pnb62 = (float) cnb62 / alignedRegionLength * 100;

    return 1;
}

int AlignmentStatistics::core(Fasta f)
{
    ncols = checkLengths(f);

    if (0 == ncols)
    {
        cout << "Sequences must all be of the same length!" << endl;
        return 0;
    }

    meanCore = 0;
    int startCore = -1;
    int endCore = -1;
    
    float nseqs = f.count();

    Sequence s = f.getSequence(0);
    int ncol = s.length();
    
    int matchcol = 0;
    int i,j;
    int matchOffset[ncol];

    for (i=0; i<ncol; i++)
    {
        char aa = s.getResidue(i).c_str()[0];
        
        if (aa == '-' || isupper(aa) )
            matchOffset[matchcol++] = i;
    }

    alignedRegionLength = matchcol;

    for (i=0; i<matchcol; i++)
    {
        int gaps = 0;
        for (j=0; j< nseqs; j++)
        {
            s = f.getSequence(j);
            if (s.getResidue(matchOffset[i]).c_str()[0] == '-')
                gaps++;
        }
        
        if (gaps / nseqs <= 0.3)
        {
            startCore = i;
            break;
        }
    }

    for (i = matchcol-1; i > startCore; i--)
    {
        int gaps = 0;
        for (j=0; j< nseqs; j++)
        {
            s = f.getSequence(j);
            if (s.getResidue(matchOffset[i]).c_str()[0] == '-')
                gaps++;
        }

        if (gaps / nseqs <= 0.3)
        {
            endCore = i+1;
            break;
        }
    }
    if (endCore == -1)
        endCore = matchcol;
    ///cout << "Start of Core: " << startCore << endl << "End of Core: " << endCore << endl;
    lengthCore = endCore - startCore;
    
    meanCoreGaps = 0.0;
    
    for (j=0; j < nseqs; j++)
    {
        s = f.getSequence(j);

        int residues = 0;
        int coreResidues = 0;
        for(i = 0; i < matchcol; i++)
        {
            if (s.getResidue(matchOffset[i]).c_str()[0] != '-')
            {
                residues++;
                if ( i >= startCore && i < endCore )
                    coreResidues++;
            }
            else
            {
                if (i >= startCore && i < endCore )
                    meanCoreGaps++;
            }
        }
//cout << "Residues: " << residues << endl
//     << "Core Residues: " << coreResidues << endl;

        meanCore += (float) coreResidues / residues;
///cout << meanCore << endl;
    }
    meanCore /= nseqs;
    meanCoreGaps /= lengthCore * nseqs;
    firstCoreRes = matchOffset[startCore];
    lastCoreRes = matchOffset[endCore-1];
}

int AlignmentStatistics::global(Fasta f)
{

    ncols = checkLengths(f);

    if (ncols == 0)
    {
        cout << "Sequences must all be of the same length!" << endl;
        return 0;
    }

    globalpi = 0.0;
    approxmeanpi = 0.0;

    Sequence s;
    Sequence t;

    float nseqs = f.count();

    s = f.getSequence(0);

    int ncol = s.length();
    int matchCol = 0;
    int char_a = (int) 'A';

    int ngaps = 0;
    float alignedPairs = 0;
    
    char aa;
    
    //loop over # positions
    for (int i = 0; i < ncol; i++)
    {
        aa = s.getResidue(i).c_str()[0];

        if (isupper(aa) || aa == '-')
        {
            matchCol++;

            int colgaps = 0;
            int k;
            int aacount[26];
            for (k = 0; k < 26; k++)
                aacount[k] = 0;
            
            for (int j = 0; j < nseqs; j++)
            {
                t = f.getSequence(j);
            
                aa = t.getResidue(i).c_str()[0];
        
                if (isupper(aa))
                    aacount[aa-char_a] ++;
                else if (aa == '-')
                {
                    ngaps++;
                    colgaps++;
                }
            }
            
            for (k = 0; k < 26; k++)
            {
                if (aacount[k] > 1)
                {
                    int ident = (aacount[k] * aacount[k] - aacount[k]) / 2;

                    //cout << ident << " ";
                    globalpi += ident;

                }
            }
            ///cout << endl;
            alignedPairs += (nseqs - colgaps) * (nseqs - colgaps -1) / 2 + colgaps * (nseqs - colgaps);
        }
    }

/*    float pctgap = ngaps / (matchCol * nseqs);
    float meanLength = matchCol * (1-pctgap);
    
    float probGap = (pctgap / nseqs / matchCol);
    probGap *= probGap;

    float approxAlignedPairs = matchCol * nseqs * (nseqs -1) / 2;
    approxAlignedPairs -= approxAlignedPairs * probGap;
    
    cout << "Number of match colummns " << matchCol << endl
        << "Number of identities " << globalpi << endl
        << "Number of gaps " << ngaps << endl
        << "Mean length of aligned region " << meanLength << endl
        << "Approx mean local percent identity " << globalpi / (meanLength * nseqs * (nseqs-1) / 2) << endl
        << "Probability of gap " << probGap << endl
        << "Approx number of aligned pairs " << approxAlignedPairs << endl
        << "Approx mean local percent identity " << globalpi / approxAlignedPairs << endl
        << "Number of aligned pairs " << alignedPairs << endl
        << "Approx mean local percent identity " << globalpi / alignedPairs << endl;
*/
    approxmeanpi = globalpi / alignedPairs;
    globalpi /= (float) (matchCol * nseqs * (nseqs-1) / 2);

}

int AlignmentStatistics::checkLengths(Fasta f)
{
    if (f.count() == 1)
    {
        return f.getSequence(0).length();
    }

    Sequence s = f.getSequence(0);

    int len = s.length();

    for (int i = 1; i < f.count(); i++)
    {
        s = f.getSequence(i);

        if (s.length() != len)
        {
            return 0;
        }
    }

    return len;
}

int AlignmentStatistics::alignedColumns()
{
    return alignedRegionLength;
}

float AlignmentStatistics::averageLength()
{
    return (float) sumSlen / (float) nseqs;
}

float AlignmentStatistics::meanAlignedToCore()
{
    return meanCore;
}

int AlignmentStatistics::coreLength()
{
    return lengthCore;
}

int AlignmentStatistics::firstCoreResidue()
{
    return firstCoreRes;
}

int AlignmentStatistics::lastCoreResidue()
{
    return lastCoreRes;
}

float AlignmentStatistics::meanCoreGapPercentage()
{
    return meanCoreGaps;
}

float AlignmentStatistics::averageAlignedRegionLength()
{
    return (float) sumSAlignedRegionLen / (float) nseqs;
}

int AlignmentStatistics::shortestLength()
{
    return shortestLen;
}

int AlignmentStatistics::shortestAlignedRegionLength()
{
    return shortestAlignedRegionLen;
}

int AlignmentStatistics::longestLength()
{
    return longestLen;
}

int AlignmentStatistics::longestAlignedRegionLength()
{
    return longestAlignedRegionLen;
}

int AlignmentStatistics::numberOfSequences()
{
    return nseqs;
}

int AlignmentStatistics::alignmentLength()
{
    return ncols;
}

string AlignmentStatistics::shortestSequenceId()
{
    return shortest;
}

string AlignmentStatistics::shortestAlignedRegionSequenceId()
{
    return shortestAlignedRegion;
}

string AlignmentStatistics::longestSequenceId()
{
    return longest;
}

string AlignmentStatistics::longestAlignedRegionSequenceId()
{
    return longestAlignedRegion;
}

float AlignmentStatistics::meanPercentIdentity()
{
    return meanpi;
}

float AlignmentStatistics::minimumPercentIdentity()
{
    return minpi;
}

string AlignmentStatistics::minimumPercentIdentityIds()
{
    return minpiids;
}

float AlignmentStatistics::maximumPercentIdentity()
{
    return maxpi;
}

float AlignmentStatistics::globalPercentIdentity()
{
    return globalpi;
}

float AlignmentStatistics::approximateMeanPercentIdentity()
{
    return approxmeanpi;
}

string AlignmentStatistics::maximumPercentIdentityIds()
{
    return maxpiids;
}

float AlignmentStatistics::percentGaps()
{
    return pgaps;
}

float AlignmentStatistics::percentColumnsNegativeBlosum62()
{
    return pnb62;
}

int AlignmentStatistics::columnsNegativeBlosum62()
{
    return cnb62;
}

int AlignmentStatistics::maximumLengthDifference()
{
    return longestLen - shortestLen;
}

int AlignmentStatistics::maximumAlignedRegionLengthDifference()
{
    return longestAlignedRegionLen - shortestAlignedRegionLen;
}

float AlignmentStatistics::percentMaximumLengthDifference()
{
    return (float) shortestLen / longestLen * 100;
}

float AlignmentStatistics::percentAlignedRegionMaximumLengthDifference()
{
    return (float) shortestAlignedRegionLen / longestAlignedRegionLen * 100;
}

void AlignmentStatistics::allowGaps(bool b)
{
    letGaps = b;
}

bool AlignmentStatistics::allowingGaps()
{
    return letGaps;
}
