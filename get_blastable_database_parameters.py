#!/usr/bin/python

import struct, glob, string

default_db = "/clusterfs/ohana/external/UniProt/current/protein"

# This parses the beginning of the index file produced by formatdb,
# as documented at
# ftp://ftp.ncbi.nlm.nih.gov/blast/documents/developer/readdb.txt

def get_blastable_database_parameters(db=default_db):
  idx_files = glob.glob("%s*.pin" % db)
  total_num_sequences = 0
  total_database_length = 0
  length_longest = 0
  for idx_file in idx_files:
    f = open(idx_file, "rb")
    ver_num = f.read(4)
    is_protein = f.read(4)
    title_len_bytes = f.read(4)
    title_len = struct.unpack(">I", title_len_bytes)
    title = f.read(title_len[0])
    dt_len_bytes = f.read(4)
    dt_len = struct.unpack(">I", dt_len_bytes)
    date_time = f.read(dt_len[0])
    num_seqs_bytes = f.read(4)
    num_seqs = struct.unpack(">L", num_seqs_bytes)
    total_num_sequences = total_num_sequences + num_seqs[0]
    db_length_bytes = f.read(8)
    db_length = struct.unpack(">L", db_length_bytes[0:4])
    total_database_length = total_database_length + db_length[0]
    len_longest_bytes = f.read(4)
    len_longest = struct.unpack(">I", len_longest_bytes)
    if len_longest[0] > length_longest:
      length_longest = len_longest[0]
  return total_num_sequences, total_database_length, length_longest

def main():
  num_sequences, db_length, length_longest = get_blastable_database_parameters()
  print "%d" % db_length

if __name__ == '__main__':
  main()

