
#include <iostream>
#include <iomanip>
using namespace std;

#include "AlignmentStatistics.h"

void usage();

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        usage();
        return(0);
    }

    bool verbose_b = true;
    string filename;
    bool bGlobal = false;
    bool bCore = false;
    
    string argv1 = argv[1];
    if (argv1 == "-h")
    {
        usage();
        return(0);
    }
    else if (argv1 == "-t" || argv1 == "-v" )
    {
        if (argc < 3)
        {
            usage();
            return(0);
        }
        else
        {
            if ( argv1 == "-t" )
            {
                verbose_b = false;
            }
            filename = argv[2];
        }
    }
    else if ( argv1 == "-g")
    {
        bGlobal = true;
        filename = argv[2];
    }
    else if ( argv1 == "-c" )
    {
        bCore = true;
        filename = argv[2];
    }
    else
    {
        filename = argv[1];
        if ( filename[0] == '-' )
        {
            cout << "Unknown option " << filename << endl;
            return( 1 );
        }
    }

    Fasta f;
    int read_err = f.readFile(filename);
    if ( read_err == 1 )
    {
       return( 1 );
    }

    AlignmentStatistics stats;

    if (bGlobal)
    {
        if (! stats.global(f) )
        {
            return 1;
        }

        cout << "GLOBALPWID:       " << stats.globalPercentIdentity() * 100 
                                                                        << endl
             << "MEANRESIDUEPWID:  " << 
                                    stats.approximateMeanPercentIdentity() * 100
                                                                        << endl;
        return 0;
    }

    if (bCore)
    {
        if (! stats.core(f) )
        {
            return 1;
        }

        cout << "MEAN_CORE_ALIGNED: " << stats.meanAlignedToCore() * 100 << endl
             << "FIRST_CORE_RES:    " << stats.firstCoreResidue()+1 << endl
             << "LAST_CORE_RES:     " << stats.lastCoreResidue()+1 << endl
             << "CORE_LENGTH:       " << stats.coreLength() << endl
             << "MEAN_CORE_GAP_PCT: " << stats.meanCoreGapPercentage() * 100 << endl
             << "ALIGNED_COLS:      " << stats.alignedColumns() << endl;
        
        return 0;
    }
    
    if (! stats.compute(f))
    {
        return 1;
    }

    if ( verbose_b )
    {
        cout 
         << setiosflags(ios::fixed)

         /* Percentages will show one decimal place (tenth of percent) */

         << setprecision(1)

         << "Input alignment filename                                   " 
         << filename
         << endl

         << "Number of sequences                                        " 
         /* Right-justify numbers. */

         << setw(6)

         << stats.numberOfSequences()
         << endl

         << "------------------------------------------------------------------"
         << endl

         << "Number of columns (alignment length)                       " 
         << setw(6)
         << stats.alignmentLength()
         << endl

         << "Average length                                             " 
         << setw(6)
         << stats.averageLength()
         << endl

         << "Shortest sequence ID                                       " 
         << stats.shortestSequenceId() 
         << endl

         << "Length of the shortest sequence                            " 
         << setw(6)
         << stats.shortestLength() 
         << endl

         << "Longest sequence ID                                        " 
         << stats.longestSequenceId() 
         << endl

         << "Length of the longest sequence                             " 
         << setw(6)
         << stats.longestLength() 
         << endl

         << "Difference in length longest-shortest                      "
         << setw(6)
         << stats.maximumLengthDifference()
         << endl

         << "Minimum overlap (fraction longest covered by shortest)     "
         << setw(6)
         << stats.percentMaximumLengthDifference() << "%"
         << endl

         << "------------------------------------------------------------------"
         << endl

         << "Number of aligned columns                                  "
         << setw(6)
         << stats.alignedColumns()
         << endl

         << "Average length (aligned columns only)                      "
         << setw(6)
         << stats.averageAlignedRegionLength()
         << endl

         << "Shortest sequence (aligned columns only) ID                " 
         << stats.shortestAlignedRegionSequenceId() 
         << endl

         << "Length of the shortest sequence (aligned columns only)     " 
         << setw(6)
         << stats.shortestAlignedRegionLength() 
         << endl

         << "Longest sequence ID (aligned columns only)                 " 
         << stats.longestAlignedRegionSequenceId() 
         << endl

         << "Length of the longest sequence (aligned columns only)      " 
         << setw(6)
         << stats.longestAlignedRegionLength() 
         << endl

         << "Difference in length longest-shortest (aligned columns)    "
         << setw(6)
         << stats.maximumAlignedRegionLengthDifference()
         << endl

         << "Minimum overlap (fraction longest covered by shortest)     "
         << setw(6)
         << stats.percentAlignedRegionMaximumLengthDifference() << "%"
         << endl

         << "------------------------------------------------------------------"
         << endl

         << "Minimum percent identity                                   " 
         << setw(6)
         << stats.minimumPercentIdentity() << "%"
         << endl

         << "Maximum percent identity                                   " 
         << setw(6)
         << stats.maximumPercentIdentity() << "%"
         << endl

         << "Mean percent identity                                      " 
         << setw(6)
         << stats.meanPercentIdentity() << "%"
         << endl

         << "Percentage of gaps over entire alignment                   " 
         << setw(6)
         << stats.percentGaps() << "%"
         << endl

         << "Percent of cols with average BLOSUM 62 scores < 0          "
         << setw(6)
         << stats.percentColumnsNegativeBlosum62() << "%"
         << endl

         << "Number of cols with average BLOSUM 62 scores < 0           "
         << setw(6)
         << stats.columnsNegativeBlosum62() 
         << endl;

    }
    else
    {
        cout << "FILE:               " << filename << endl
             << "SEQS:               " << stats.numberOfSequences() << endl
             << "COLS:               " << stats.alignmentLength() << endl
             << "AVERAGE_LENGTH:     " << stats.averageLength() << endl
             << "SHORTEST:           " << stats.shortestSequenceId() << endl
             << "SHORTEST_LENGTH:    " << stats.shortestLength() << endl
             << "LONGEST:            " << stats.longestSequenceId() << endl
             << "LONGEST_LENGTH:     " << stats.longestLength() << endl
             << "MAXLENDIFF:         " << stats.maximumLengthDifference() 
                                                                         << endl
             << "PMINOVER:           " << stats.percentMaximumLengthDifference() 
                                                                         << endl

             << "ALIGNED_COLS:       " << stats.alignedColumns() <<endl
             << "AVERAGE_LENGTH_AC:  " << stats.averageAlignedRegionLength() 
                                                                         << endl
             << "SHORTEST_AC:        " 
                                      << stats.shortestAlignedRegionSequenceId()
                                                                         << endl
             << "SHORTEST_LENGTH_AC: " << stats.shortestAlignedRegionLength() 
                                                                         << endl
             << "LONGEST_AC:         " << stats.longestAlignedRegionSequenceId()
                                                                         << endl
             << "LONGEST_LENGTH_AC:  " << stats.longestAlignedRegionLength() 
                                                                         << endl
             << "MAXLENDIFF_AC:      " 
                                 << stats.maximumAlignedRegionLengthDifference()
                                                                         << endl
             << "PMINOVER_AC:        " 
                          << stats.percentAlignedRegionMaximumLengthDifference()
                                                                         << endl

             << "MINPID:             " << stats.minimumPercentIdentity() << endl
             << "MAXPID:             " << stats.maximumPercentIdentity() << endl
             << "MEANPID:            " << stats.meanPercentIdentity() << endl
             << "PGAPS:              " << stats.percentGaps() << endl
             << "PCNAB62:            " << stats.percentColumnsNegativeBlosum62()
                                                                         << endl
             << "NCNAB62:            " << stats.columnsNegativeBlosum62() 
                                                                        << endl;
    }

    return(0);
}

void usage()
{
    cout << endl
         << "Usage: astats [ -t or -g or -c ] <msafile>" << endl
         << endl
         << "Options:" << endl
         << "  -t  Terse output." << endl
         << "  -g  Global percent ID quick calculation." << endl
         << "  -c  Determine core aligned region of alignment." << endl
         << endl
         << "Berkeley Phylogenomics Group - http://phylogenomics.berkeley.edu"
                                                                        << endl
         << "Copyright(c) 2005 Regents of the University of California"
                                                                        << endl
         << endl;
}


