#!/usr/bin/python
import os, sys, string
from optparse import OptionParser

def main():
  usage = """%s [options] <organism1> <organism2> <organism3> ...
  
  Here <organism1> is the five character SwissProt identifier 
  for the organism, such as CHICK for chicken or DROME for fruitfly.""" \
  % sys.argv[0]
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-e", "--max_e_value", dest="evalue_cutoff",
                        default="1e-3", 
                        help="Maximum e-value for which to include hits")
  (options, args) = opt_parser.parse_args()
  if len(args) == 0:
    opt_parser.error("No organisms were specified")
  try:
    evalue_cutoff = float(options.evalue_cutoff)
  except ValueError:
    opt_parser.error("--max_e_value must be a number")

  organisms = []
  for arg in args:
    if len(arg) == 5 and arg.isalpha():
      organisms.append(string.upper(arg))
    else:
      opt_parser.error("%s is not a 5-character alphabetical organism id" % arg)

  all_hits = {}
  hits_by_organism = {}
  for organism in organisms:
    hits_by_organism[organism] = {}

  patt = '|'.join(['_%s' % organism for organism in organisms])
  # The following would retrieve the identifiers of all the PSI-BLAST hits
  # cmd = 'grep -E "%s" pb | cut -f 2 | cut -d "|" -f 3' % patt
  cmd = 'grep -E "%s" pb' % patt
  fp_grep_cmd = 'grep -E "%s" final.a2m' % patt
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
      f = os.popen(fp_grep_cmd)
      fp_hit_lines = f.readlines()
      f.close()
      fp_hit_set = set()
      for fp_hit_line in fp_hit_lines:
        fp_hit_id = fp_hit_line.rstrip().strip('>').split()[0]
        if fp_hit_id.find('|') >= 0:
          if fp_hit_id[0:3] == 'sp|' or fp_hit_id[0:4] == 'lcl|':
            fp_hit_id = fp_hit_id.split('|')[2]
          else:
            fp_hit_id = fp_hit_id.split('|')[1]
        fp_hit_set.add(fp_hit_id)
      for hit_line in hit_lines:
        hit_line = hit_line.rstrip()
        fields = hit_line.split()
        evalue = float(fields[10])
        if evalue > evalue_cutoff:
          continue
        hit_id = fields[1].split('|')[2]
        is_fp_hit = hit_id in fp_hit_set
        organism = hit_id.split('_')[1]
        if hit_id not in hits_by_organism[organism]:
          hits_by_organism[organism][hit_id] = (evalue, seed_id, is_fp_hit)
          all_hits[hit_id] = (evalue, seed_id, is_fp_hit)
        elif hits_by_organism[organism][hit_id][0] > evalue:
          hits_by_organism[organism][hit_id] = (evalue, seed_id, is_fp_hit)
          all_hits[hit_id] = (evalue, seed_id, is_fp_hit)
      os.chdir('..')
      found_seed = True
  print "Organisms:", organisms
  hits_with_evalues = [(all_hits[hit], hit) for hit in all_hits.keys()]
  hits_with_evalues.sort()
  for (evalue, seed_id, is_fp_hit), hit_id in hits_with_evalues:
    if is_fp_hit:
      print "%s\t%g\t%s\t*" % (hit_id, evalue, seed_id)
    else:
      print "%s\t%g\t%s" % (hit_id, evalue, seed_id)


if __name__ == '__main__':
  main()
