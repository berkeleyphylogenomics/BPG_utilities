#!/usr/bin/env python

from optparse import OptionParser
import os
from stat import ST_SIZE
import sys


"""Evenly divide the UniProt file on the UniProt record boundaries.

Over time, the UniProt.dat data file is exponentially increasing in size
as scientists discover more information. To make the job of inserting
this data into the database more managable, we need to evenly divide the
file into discrete chunks.  Because a single UniProt record spans
multiple rows in the file, care needs to be taken to split the file
along the UniProt record boundaries.

A record begins wih the ID line and ends with the sentinal value of
'//\n'. For example:

ID   002R_IIV3               Reviewed;         458 AA.
AC   Q197F8;
DT   16-JUN-2009, integrated into UniProtKB/Swiss-Prot.

[snip]

     STILHKRDTD WVENNPLKTP AQVEMYKFLL RISQLNRDGT GYESDSDPEN EHFDDESFSS
     GEEDSSDEDD PTWAPDSDDS DWETETEEEP SVAARILEKG KLTITNLMKS LGFKPKPKKI
     QSIDRYFCSL DSNYNSEDED FEYDSDSEDD DSDSEDDC
//

So, divisions should be immediate after the sentinel value '//\n'.

The file does not need to be physically divided into chunks, as the
load_uniprot.py program can understand where to start in a file and
where to stop reading. For example, if a record started at 30930627407
bytes into the file, load_uniprot.py can start reading there.
Additionally, if the UniProt record ended at byte offset 31246245389
within the file, the load_uniprot.py program can stop reading at that
point.

Thus, we need to tell load_uniprot.py where to start reading and where
to stop reading within the file. In this example mentioned above, the
command would look like this:

load_uniprot.py 30930627407_31246245389

This script takes the number of divisions requested, opens the
uniprot.dat file to import, and calculates the offset boundaries. 

Example outputs from the program follow:

315621241_631238004
631238005_946857079
946857080_1262475559
1262475560_1578094527
1578094528_1893714137
1893714138_2209330787
2209330788_2524951107
2524951108_2840569614


Each of these values can be given to the load_uniprot.py script, so many
instances of the program can be running simultaneously.
"""


def main(options, args):

  uniprot_dat_path =\
               "/clusterfs/ohana/external/UniProt/to_import/uniprot.dat"
  uniprot_dat_size = os.stat(uniprot_dat_path)[ST_SIZE]

  # Save the last shard for later
  target_shard_size = uniprot_dat_size / (num_shards - 1)

  f = open(uniprot_dat_path)
  shard_start_pos = f.tell()

  for shard_num in xrange(num_shards-1):
    # We divide the into the exact number of shards. But, many records
    # would be divided in half. Instead, extend the end of the recod to
    # be the end of the UniProt record.

    # Seek into position
    f.seek(shard_num*target_shard_size)

    # Advance until the end of the UniProt record that has been started
    line = f.readline()
    while line != '//\n':
        line = f.readline()

    # Save new start position and print results
    sys.stdout.write("%d_%d\n" % (shard_start_pos, f.tell()-1))
    shard_start_pos = f.tell()

  # The last shard takes the remaining portion to the end of the file
  sys.stdout.write("%d_%d\n" % (shard_start_pos, uniprot_dat_size))



if __name__ == '__main__':

  parser = OptionParser(usage='%prog')
  parser.add_option(
      "-n", "--num_shards", dest="num_shards",
      default="100", 
      help="number of shards into which to split uniprot.dat file")
  (options, args) = parser.parse_args()

  try:
    num_shards = int(options.num_shards)
  except ValueError:
    parser.error("num_shards must be an integer")

  if num_shards <= 0:
    parser.error("num_shards must be a positive number")

  main(options, args)
