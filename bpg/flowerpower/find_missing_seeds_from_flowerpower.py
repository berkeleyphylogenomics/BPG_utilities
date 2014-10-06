#!/usr/bin/python
import os, sys, string
from optparse import OptionParser
from Bio import Seq, SeqIO

def main():
  usage = """"%s fasta_file [organism]

  fasta_file should contain all the sequences that were originally clustered
  into seeds for FlowerPower, in FASTA format.
  %s will check that each of these sequences is present in the set of 
  homologues retrieved by FlowerPower for some seed.
  organism, if provided, should be the five-character SwissProt/TrEMBL 
  identifier of the organism whose genome was clustered.
  """ % (sys.argv[0], sys.argv[0])
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) < 1:
    opt_parser.error("Please provide the fasta_file argument." % sys.argv[0])
  organism=""
  if len(args) >= 2:
    if len(args[1]) == 5 and args[1].isalpha():
      organism=string.upper(args[1])
    else:
      opt_parser.error("%s is not a 5-character alphabetical organism id" \
                        % args[1])

  if organism == "":
    patt = "^>"
  else:
    patt = "^>[^ ]*_%s" % organism
  cmd = "grep -E \"%s\" final.a2m | cut -d ' ' -f 1" % patt 
  universe_handle = open(args[0], "rU")
  universe = SeqIO.to_dict(SeqIO.parse(universe_handle, "fasta"))
  universe_handle.close()
  f = open("clusters", "rU")
  lines = f.readlines()
  f.close()
  cluster_num = 0
  found_seed = False
  for line in lines:
    line = line.rstrip()
    if line[0:9] == 'Cluster #':
      cluster_num = int(line[9:-1])
      found_seed = False
    elif not found_seed:
      if line.find('|') >= 0:
        seed_id = line.split('|')[1]
      else:
        seed_id = line
      # Replace / in the header by %2f
      seed_id = seed_id.replace('/','%2f')
      # Replace ( in the header by %28
      seed_id = seed_id.replace('(','%28')
      # Replace ) in the header by %29
      seed_id = seed_id.replace(')','%29')
      # Replace - in the header by %2d
      seed_id = seed_id.replace('-','%2d')
      # Replace * in the header by %2a
      seed_id = seed_id.replace('*','%2a')
      os.chdir(seed_id)
      f = os.popen(cmd)
      hit_lines = f.readlines()
      f.close()
      for hit_line in hit_lines:
        hit_id = hit_line.rstrip().strip('>').split()[0]
        print "Found hit %s" % hit_id
        if hit_id in universe:
          print "hit was in universe, removing"
          del universe[hit_id]
      os.chdir('..')
      found_seed = True
  if len(universe.keys()) > 0:
    for hit_id in universe.keys():
      print hit_id
  else:
    print "All sequences were retrieved!"

if __name__ == '__main__':
  main()
