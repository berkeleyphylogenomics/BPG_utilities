#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

void usage(char *program_name) {
  printf( "Usage: %s [--assume_valid_a2m] <alignment_file_path>\n\n",
            program_name );

  printf( "This program takes as input a protein MSA (multiple sequence " );
  printf( "file in a2m (see \n" );
  printf( "http://compbio.soe.ucsc.edu/a2m-desc.html) format.  It requires " );
  printf( "the alignment\n" );
  printf( "to be prettyaligned (i.e., with dots inserted so all aligned " );
  printf( "sequences are the\n" );
  printf( "same length).\n\n" );
  printf( "The program prints comma-separated rows with the column index " );
  printf( "(zero-based) of\n" );
  printf( "each aligned column with respect to the overall prettyaligned " );
  printf( "alignment; the\n" );
  printf( "most conserved residue in that column; and the BLOSUM62 " );
  printf( "sum-of-pairs score\n" );
  printf( "for the column.\n\n" );

  printf( "If you know the input is a valid a2m file (e.g., because another " );
  printf( "program has \n" );
  printf( "already validated it), you can issue the --assume_valid_a2m " );
  printf( "option to \n" );
  printf( "skip validity checking.\n\n" );
  
}

void invalid_a2m(char *alignment_file_path, char *msg) {
  fprintf( stderr, "%s is not a valid a2m file.  %s Exiting...\n",
          alignment_file_path, msg );
  exit(1);
}

void malloc_failed(char *variable_name) {
  fprintf( stderr, "Unable to allocate memory for %s. Exiting...\n",
          variable_name );
  exit(1);
}

int main(int argc, char **argv) {
  char *alignment_file_path;
  int check_a2m = 1;

  switch (argc) {
    case 2:
      if (*(argv[1]) == '-') {
        usage(argv[0]);
        exit(0);
      } else {
        alignment_file_path = argv[1];
      }
      break;
    case 3:
      if (strcmp(argv[1], "--assume_valid_a2m") || *(argv[2]) == '-') {
        usage(argv[0]);
        exit(0);
      } else {
        check_a2m = 0;
        alignment_file_path = argv[2];
      }
      break;
    default:
      usage(argv[0]);
      exit(0);
      break;
  };

  struct stat alignment_file_stat;

  if (stat(alignment_file_path, &alignment_file_stat) != 0) {
    fprintf( stderr, "Error: %s\n", strerror(errno) );
    fprintf( stderr, "Unable to stat %s; exiting...\n", alignment_file_path );
    exit(errno);
  }

  size_t file_size = alignment_file_stat.st_size;

  char *alignment_bytes = (char *)malloc(file_size + 1);

  if (!alignment_bytes) {
    malloc_failed("alignment_bytes");
  }

  FILE *alignment_f = fopen(alignment_file_path, "r");

  if (!alignment_f) {
    fprintf( stderr, "Error: %s\n", strerror(errno) );
    fprintf( stderr, "Unable to open %s; exiting...\n", alignment_file_path );
    exit(errno);
  }

  size_t bytes_read = fread(alignment_bytes, 1, file_size, alignment_f);

  fclose(alignment_f);

  if (check_a2m && *alignment_bytes != '>') {
    invalid_a2m( alignment_file_path, "It doesn't start with '>'." );
  }

  alignment_bytes[bytes_read] = '\0';

  int num_sequences = 0;

  char *cur_ptr;

  char *first_seq_ptr = NULL;

  for (cur_ptr = alignment_bytes; *cur_ptr; ++cur_ptr) {
    if (*cur_ptr == '>') {
      ++num_sequences;
    } else if (!first_seq_ptr && *cur_ptr == '\n') {
      first_seq_ptr = ++cur_ptr;
    }
  }


  if (num_sequences == 0 || !first_seq_ptr) {
    fprintf( stderr, "No sequences found.  Exiting...\n" );
    exit(1);
  } 

  int j = 0;

  /* Declare the BLOSUM62 matrix.  This particular version of the matrix is from
   * http://icb.med.cornell.edu/education/courses/introbio/BLOSUM62, which
   * includes entries for the letters B, J, Z, and X.  For the letters O and U,
   * which do not code for amino acids (even ambiguously), we use the value -4.
   * We use doubles since we will need them anyway for the sum-of-pairs scores.
   */
  double blosum62[26][26] = { 
   /*     A,    B,    C,    D,    E,    F,    G,    H,    I,
    *     J,    K,    L,    M,    N,    O,    P,    Q,    R,
    *     S,    T,    U,    V,    W,    X,    Y,    Z
    */
    /* A */
    {   4.0, -2.0,  0.0, -2.0, -1.0, -2.0,  0.0, -2.0, -1.0,
       -1.0, -1.0, -1.0, -1.0, -2.0, -4.0, -1.0, -1.0, -1.0,
        1.0,  0.0, -4.0,  0.0, -3.0, -1.0, -2.0, -1.0 },
    /* B */
    {  -2.0,  4.0, -3.0,  4.0,  1.0, -3.0, -1.0,  0.0, -3.0,
       -3.0,  0.0, -4.0, -3.0,  4.0, -4.0, -2.0,  0.0, -1.0,
        0.0, -1.0, -4.0, -3.0, -4.0, -1.0, -3.0,  0.0 },
    /* C */
    {   0.0, -3.0,  9.0, -3.0, -4.0, -2.0, -3.0, -3.0, -1.0,
       -1.0, -3.0, -1.0, -1.0, -3.0, -4.0, -3.0, -3.0, -3.0,
       -1.0, -1.0, -4.0, -1.0, -2.0, -1.0, -2.0, -3.0 },
    /* D */
    {  -2.0,  4.0, -3.0,  6.0,  2.0, -3.0, -1.0, -1.0, -3.0,
       -3.0, -1.0, -4.0, -3.0,  1.0, -4.0, -1.0,  0.0, -2.0,
        0.0, -1.0, -4.0, -3.0, -4.0, -1.0, -3.0,  1.0 },
    /* E */
    {  -1.0,  1.0, -4.0,  2.0,  5.0, -3.0, -2.0,  0.0, -3.0,
       -3.0,  1.0, -3.0, -2.0,  0.0, -4.0, -1.0,  2.0,  0.0,
        0.0, -1.0, -4.0, -2.0, -3.0, -1.0, -2.0,  4.0 },
    /* F */
    {  -2.0, -3.0, -2.0, -3.0, -3.0,  6.0, -3.0, -1.0,  0.0,
        0.0, -3.0,  0.0,  0.0, -3.0, -4.0, -4.0, -3.0, -3.0,
       -2.0, -2.0, -4.0, -1.0,  1.0, -1.0,  3.0, -3.0 },
    /* G */
    {   0.0, -1.0, -3.0, -1.0, -2.0, -3.0,  6.0, -2.0, -4.0,
       -4.0, -2.0, -4.0, -3.0,  0.0, -4.0, -2.0, -2.0, -2.0,
        0.0, -2.0, -4.0, -3.0, -2.0, -1.0, -3.0, -2.0 },
    /* H */
    {  -2.0,  0.0, -3.0, -1.0,  0.0, -1.0, -2.0,  8.0, -3.0,
       -3.0, -1.0, -3.0, -2.0,  1.0, -4.0, -2.0,  0.0,  0.0,
       -1.0, -2.0, -4.0, -3.0, -2.0, -1.0,  2.0,  0.0 },
    /* I */
    {  -1.0, -3.0, -1.0, -3.0, -3.0,  0.0, -4.0, -3.0,  4.0,
        3.0, -3.0,  2.0,  1.0, -3.0, -4.0, -3.0, -3.0, -3.0,
       -2.0, -1.0, -4.0,  3.0, -3.0, -1.0, -1.0, -3.0 },
    /* J */
    {  -1.0, -3.0, -1.0, -3.0, -3.0,  0.0, -4.0, -3.0,  3.0,
        3.0, -3.0,  3.0,  2.0, -3.0, -4.0, -3.0, -2.0, -2.0,
       -2.0, -1.0, -4.0,  2.0, -2.0, -1.0, -1.0, -3.0 },
    /* K */
    {  -1.0,  0.0, -3.0, -1.0,  1.0, -3.0, -2.0, -1.0, -3.0,
       -3.0,  5.0, -2.0, -1.0,  0.0, -4.0, -1.0,  1.0,  2.0,
        0.0, -1.0, -4.0, -2.0, -3.0, -1.0, -2.0,  1.0 },
    /* L */
    {  -1.0, -4.0, -1.0, -4.0, -3.0,  0.0, -4.0, -3.0,  2.0,
        3.0, -2.0,  4.0,  2.0, -3.0, -4.0, -3.0, -2.0, -2.0,
       -2.0, -1.0, -4.0,  1.0, -2.0, -1.0, -1.0, -3.0 },
    /* M */
    {  -1.0, -3.0, -1.0, -3.0, -2.0,  0.0, -3.0, -2.0,  1.0,
        2.0, -1.0,  2.0,  5.0, -2.0, -4.0, -2.0,  0.0, -1.0,
       -1.0, -1.0, -4.0,  1.0, -1.0, -1.0, -1.0, -1.0 },
   /*     A,    B,    C,    D,    E,    F,    G,    H,    I,
    *     J,    K,    L,    M,    N,    O,    P,    Q,    R,
    *     S,    T,    U,    V,    W,    X,    Y,    Z
    */
    /* N */
    {  -2.0,  4.0, -3.0,  1.0,  0.0, -3.0,  0.0,  1.0, -3.0,
       -3.0,  0.0, -3.0, -2.0,  6.0, -4.0, -2.0,  0.0,  0.0,
        1.0,  0.0, -4.0, -3.0, -4.0, -1.0, -2.0,  0.0 },
    /* O */
    {  -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,
       -4.0, -4.0, -4.0, -4.0, -4.0,  1.0, -4.0, -4.0, -4.0,
       -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0 },
    /* P */
    {  -1.0, -2.0, -3.0, -1.0, -1.0, -4.0, -2.0, -2.0, -3.0,
       -3.0, -1.0, -3.0, -2.0, -2.0, -4.0,  7.0, -1.0, -2.0,
       -1.0, -1.0, -4.0, -2.0, -4.0, -1.0, -3.0, -1.0 },
    /* Q */
    {  -1.0,  0.0, -3.0,  0.0,  2.0, -3.0, -2.0,  0.0, -3.0,
       -2.0,  1.0, -2.0,  0.0,  0.0, -4.0, -1.0,  5.0,  1.0,
        0.0, -1.0, -4.0, -2.0, -2.0, -1.0, -1.0,  4.0 },
    /* R */
    {  -1.0, -1.0, -3.0, -2.0,  0.0, -3.0, -2.0,  0.0, -3.0,
       -2.0,  2.0, -2.0, -1.0,  0.0, -4.0, -2.0,  1.0,  5.0,
       -1.0, -1.0, -4.0, -3.0, -3.0, -1.0, -2.0,  0.0 },
    /* S */
    {   1.0,  0.0, -1.0,  0.0,  0.0, -2.0,  0.0, -1.0, -2.0,
       -2.0,  0.0, -2.0, -1.0,  1.0, -4.0, -1.0,  0.0, -1.0,
        4.0,  1.0, -4.0, -2.0, -3.0, -1.0, -2.0,  0.0 },
    /* T */
    {   0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -2.0, -1.0,
       -1.0, -1.0, -1.0, -1.0,  0.0, -4.0, -1.0, -1.0, -1.0,
        1.0,  5.0, -4.0,  0.0, -2.0, -1.0, -2.0, -1.0 },
    /* U */
    {  -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,
       -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,
       -4.0, -4.0,  1.0, -4.0, -4.0, -4.0, -4.0, -4.0 },
    /* V */
    {   0.0, -3.0, -1.0, -3.0, -2.0, -1.0, -3.0, -3.0,  3.0,
        2.0, -2.0,  1.0,  1.0, -3.0, -4.0, -2.0, -2.0, -3.0,
       -2.0,  0.0, -4.0,  4.0, -3.0, -1.0, -1.0, -2.0 },
    /* W */
    {  -3.0, -4.0, -2.0, -4.0, -3.0,  1.0, -2.0, -2.0, -3.0,
       -2.0, -3.0, -2.0, -1.0, -4.0, -4.0, -4.0, -2.0, -3.0,
       -3.0, -2.0, -4.0, -3.0, 11.0, -1.0,  2.0, -2.0 },
    /* X */
    {  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
       -1.0, -1.0, -1.0, -1.0, -1.0, -4.0, -1.0, -1.0, -1.0,
       -1.0, -1.0, -4.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
    /* Y */
    {  -2.0, -3.0, -2.0, -3.0, -2.0,  3.0, -3.0,  2.0, -1.0,
       -1.0, -2.0, -1.0, -1.0, -2.0, -4.0, -3.0, -1.0, -2.0,
       -2.0, -2.0, -4.0, -1.0,  2.0, -1.0,  7.0, -2.0 },
    /* Z */
    {  -1.0,  0.0, -3.0,  1.0,  4.0, -3.0, -2.0,  0.0, -3.0,
       -3.0,  1.0, -3.0, -1.0,  0.0, -4.0, -1.0,  4.0,  0.0,
        0.0, -1.0, -4.0, -2.0, -2.0, -1.0, -2.0,  4.0 }
  };
    
  int n, m;

  for (n = 0; n < 26; ++n) {
    for (m = 0; m < n; ++m) {
      if (blosum62[n][m] != blosum62[m][n]) {
        printf( "Symmetry violation: n = %d (%c), m = %d (%c)\n",
                n, 'A' + n, m, 'A' + m );
      }
    }
  }
  
  if (num_sequences == 1) {
    /* Every column is trivially conserved. */
    printf("ColumnIndex,ConservedResidue,Blosum62ConservationScore\n");
    for (cur_ptr = first_seq_ptr; *cur_ptr; ++cur_ptr) {
      if (*cur_ptr == '-' || *cur_ptr >= 'A' && *cur_ptr <= 'Z') {
        ++j;
        if (*cur_ptr != '-') {
          printf( "%d,%c,%g\n", j, *cur_ptr, 
                blosum62[*cur_ptr-'A'][*cur_ptr-'A'] );
        }
      } else if (check_a2m && 
                *cur_ptr != '\t' && *cur_ptr != ' ' && *cur_ptr != '\n' &&
                *cur_ptr != '.' && (*cur_ptr < 'a' || *cur_ptr > 'z') ) {
        invalid_a2m( alignment_file_path,
          "Its sequence portion contains invalid character(s)." );
      }
    }
    exit(0);
  }

  char **alignment_rows = (char **)malloc(num_sequences * sizeof(char *));

  if (!alignment_rows) {
    malloc_failed("alignment_rows");
  }

  int seeking_header = 1;
  int num_to_shift = 0;
  int total_length;
  int i = -1, col_index = 0;
  /* We could save space with a set implementation as a bit array, but let's
   * keep things simple. */
  unsigned char *is_column_aligned;
  int *aligned_column_indices;
  int num_aligned_columns = 0;

  /* After this loop, alignment_rows will have null-terminated strings with the
   * aligned sequences, with whitespace stripped.  These strings are all part of
   * the memory block alignment_bytes; whitespace is stripped in place.  The
   * array 
   */

  for (cur_ptr = alignment_bytes; *cur_ptr; ++cur_ptr) {
    if (seeking_header) {
      /* We're within the sequence part of the FASTA record */
      if (*cur_ptr == '>') {
        /* We're starting a new FASTA record. */
        *(cur_ptr - num_to_shift)= '\0';  /* Null-terminate the previous one. */
        seeking_header = 0;
        if (i == 0) {
          /* We have just finished traversing the first FASTA record. */
          total_length = j;
          num_aligned_columns = col_index;
          aligned_column_indices = (int*)malloc(col_index * sizeof(int));
          if (!aligned_column_indices) {
            malloc_failed("aligned_column_indices");
          }
          if (check_a2m) {
            is_column_aligned = (unsigned char*)malloc(total_length *
                                                        sizeof(unsigned char));
            if (!is_column_aligned) {
              malloc_failed("is_column_aligned");
            }
          }
        } else if (check_a2m && i >= 1) {
          if (j != total_length) {
            invalid_a2m( alignment_file_path,
              "The sequences are of different lengths." );
          } 
          if (col_index != num_aligned_columns) {
            invalid_a2m( alignment_file_path,
              "The columns don't line up." );
          }
        }
        ++i;
        num_to_shift = 0;
        j = 0;
        col_index = 0;
      } else if (*cur_ptr == '\n' || *cur_ptr == '\t' || *cur_ptr == ' ') {
        /* Add to the number of whitespace characters we will close up. */
        ++num_to_shift;
      } else { 
        if (num_to_shift > 0) {
          /* Strip previous whitespace by copying in place. */
          *(cur_ptr - num_to_shift) = *cur_ptr;
        }
        if (check_a2m && i > 0 && j > total_length) {
          invalid_a2m(alignment_file_path, 
            "The sequences are of different lengths.");
        }
        if (*cur_ptr == '-' || *cur_ptr >= 'A' && *cur_ptr <= 'Z') {
          if (i == 1) {
            aligned_column_indices[col_index] = j;
          }
          if (check_a2m && i >= 1) {
            if (col_index > num_aligned_columns) {
              invalid_a2m(alignment_file_path,
                "The columns don't line up.");
            }
            if (i == 1) {
              is_column_aligned[j] = 1;
            } else {
              if (!is_column_aligned[j]) {
                invalid_a2m(alignment_file_path,
                  "The columns don't line up.");
              }
            }
          }
          ++col_index;
        } else if (check_a2m) {
          if (*cur_ptr != '.' && (*cur_ptr < 'a' || *cur_ptr > 'z')) {
            invalid_a2m(alignment_file_path,
              "Its sequence portion contains invalid character(s).");
          }
          if (i == 1) {
            is_column_aligned[j] = 0;
          } else if (i > 1) {
            if (is_column_aligned[j]) {
              invalid_a2m(alignment_file_path,
                "The columns don't line up.");
            }
          }
        }
        ++j;
      }
    } else {
      /* We're within the header part of the FASTA record */
      if (*cur_ptr == '\n') {
        seeking_header = 1;
        /* We had marked the end of the file with a '\0'. */
        if (cur_ptr[1]) {
          alignment_rows[i] = cur_ptr + 1;
        }
        j = 0;
        col_index = 0;
        num_to_shift = 0;
      }
    }
  }

  if (check_a2m) {
    /* Check that the aligned characters of the i=0 row are consistent. */
    for (cur_ptr = alignment_rows[0], j = 0; *cur_ptr; ++cur_ptr, ++j) {
      if (*cur_ptr == '-' || *cur_ptr >= 'A' && *cur_ptr <= 'Z') {
        if (!is_column_aligned[j]) {
          invalid_a2m(alignment_file_path,
                  "The columns don't line up.");
        }
      } else {
        /* We already checked this is a lowercase character or dot in the
         * previous loop.
         */
        if (is_column_aligned[j]) {
          invalid_a2m(alignment_file_path,
              "The columns don't line up.");
        }
      }
    }
  }

  printf("ColumnIndex,ConservedResidue,Blosum62ConservationScore\n");
  unsigned int frequency_of_residue[26], highest_frequency;
  int residue_index;
  char residue, residue0, residue1, most_frequent_residue;
  int num_pairs, i0, i1;
  double sum_of_scores;
  for (col_index = 0; col_index < num_aligned_columns; ++col_index) {
    j = aligned_column_indices[col_index];
    for (residue_index = 0; residue_index < 26; ++residue_index) {
      frequency_of_residue[residue_index] = 0;
    }
    highest_frequency = 0;
    most_frequent_residue = '\0';
    for (i = 0; i < num_sequences; ++i) {
      residue = alignment_rows[i][j];
      if (residue == '-') {
        continue;
      }
      residue_index = residue - 'A';
      ++frequency_of_residue[residue_index];
      if ( frequency_of_residue[residue_index] > highest_frequency ) {
        highest_frequency = frequency_of_residue[residue_index];
        most_frequent_residue = residue;
      }
    }
    num_pairs = 0;
    sum_of_scores = 0.0;
    for (i0 = 0; i0 < num_sequences; ++i0) {
      residue0 = alignment_rows[i0][j];
      if (residue0 != '-') {
        for (i1 = 0; i1 < i0; ++i1) {
          residue1 = alignment_rows[i1][j];
          if (residue1 != '-') {
            sum_of_scores += blosum62[residue0-'A'][residue1-'A'];
            ++num_pairs;
          }
        }
      }
    }
    if (num_pairs > 0) {
      printf("%d,%c,%0.6f\n", j, most_frequent_residue, 
              sum_of_scores / num_pairs );
    } else {
      fprintf( stderr, "Gappy %dth aligned column %d\n", col_index, j );
    }
  }

  if (check_a2m) {
    free(is_column_aligned);
  }
  free(aligned_column_indices);
  free(alignment_rows);
  free(alignment_bytes);
}
