#!/usr/bin/python
import os, sys, string
from optparse import OptionParser
from Bio import SeqIO
import cPickle

def main():
  usage = """%s [options] fasta_file <organism1> <organism2> <organism3> ...
  
  Here id_file contains the ids of the seeds for FlowerPower.
  Here <organism1> is the five character SwissProt identifier 
  for the organism, such as CHICK for chicken or DROME for fruitfly.""" \
  % sys.argv[0]
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-e", "--max_e_value", dest="evalue_cutoff",
                        default="1e-3", 
                        help="Maximum e-value for which to include hits")
  (options, args) = opt_parser.parse_args()
  if len(args) < 2:
    opt_parser.error("No organisms were specified")
  if not os.path.exists(args[0]):
    opt_parser.error("Couldn't find %s" % args[0])
  else:
    id_file = args[0]
  try:
    evalue_cutoff = float(options.evalue_cutoff)
  except ValueError:
    opt_parser.error("--max_e_value must be a number")

  organisms = []
  for arg in args[1:]:
    arg = string.upper(arg)
    if len(arg) == 5 and arg.isalpha() or \
        len(arg) == 3 and arg in ['RAT', 'PIG', 'PEA']:
      organisms.append(arg)
    else:
      opt_parser.error("%s is not a 5-character alphabetical organism id" % arg)

  patt = '|'.join(['_%s' % organism for organism in organisms])
  # The following would retrieve the identifiers of all the PSI-BLAST hits
  # cmd = 'grep -E "%s" pb | cut -f 2 | cut -d "|" -f 3' % patt
  cmd = 'grep -E "%s" pb' % patt
  fp_grep_cmd = 'grep -e "^>" final.a2m'
  sp_grep_cmd = "grep -e 'sp|' pb_sprot"
  cur_dir = os.getcwd()
  fp_hits = {}
  f = open(id_file, "rU")
  seed_lines = f.readlines()
  f.close()
  seed_num = 0
  for line in seed_lines:
    seed_num = seed_num + 1
    line = line.rstrip()
    if line.find('|') >= 0:
      seed_id = line.split('|')[1]
    else:
      seed_id = line
    if not os.path.exists(seed_id):
      print "Currently in directory %s" % os.getcwd()
      print "Can't find %s" % seed_id
    os.chdir(seed_id)
    if not os.path.exists('pb_sprot'):
      print "pb_sprot missing from %s" % seed_id
    fp_f = os.popen(fp_grep_cmd)
    fp_hit_lines = fp_f.readlines()
    fp_f.close()
    for fp_hit_line in fp_hit_lines:
      fp_hit_id = fp_hit_line.rstrip().strip('>').split()[0]
      if fp_hit_id.find('|') >= 0:
        if fp_hit_id[0:3] == 'sp|' or fp_hit_id[0:4] == 'lcl|':
          fp_hit_id = fp_hit_id.split('|')[2]
        else:
          fp_hit_id = fp_hit_id.split('|')[1]
      if fp_hit_id in fp_hits and fp_hits[fp_hit_id] != None:
        fp_hits[fp_hit_id] = fp_hits[fp_hit_id].append(seed_id)
      else:
        fp_hits[fp_hit_id] = [seed_id]
    os.chdir(cur_dir)
  print "%d total FlowerPower hits" % len(fp_hits.keys())
  fp_pickle_f = open("FlowerPowerHits.pkl", "w")
  cPickle.dump(fp_hits, fp_pickle_f)
  fp_pickle_f.close()
  sp_hit_set = set()
  seed_num = 0
  for line in seed_lines:
    seed_num = seed_num + 1
    line = line.rstrip()
    if line.find('|') >= 0:
      seed_id = line.split('|')[1]
    else:
      seed_id = line
    if not os.path.exists(seed_id):
      print "Currently in directory %s" % os.getcwd()
      print "Can't find %s" % seed_id
    os.chdir(seed_id)
    hit_f = os.popen(cmd)
    hit_lines = hit_f.readlines()
    hit_f.close()
    sp_f = os.popen(sp_grep_cmd)
    sp_hit_lines = sp_f.readlines()
    sp_f.close()
    best_sp_hit = None
    best_sp_evalue = 100
    best_sp_hit_is_in_fp = False
    for sp_hit_line in sp_hit_lines:
      sp_hit_line = sp_hit_line.rstrip()
      fields = sp_hit_line.split()
      evalue = float(fields[10])
      sp_hit_id = fields[1].split('|')[2]
      sp_hit_set.add(sp_hit_id)
      is_fp_hit = sp_hit_id in fp_hits and fp_hits[sp_hit_id] != None
      if is_fp_hit:
        if not best_sp_hit_is_in_fp or evalue < best_sp_evalue:
          best_sp_hit_is_in_fp = True
          best_sp_hit = sp_hit_id
          best_sp_evalue = evalue
      else:
        if evalue > evalue_cutoff or best_sp_hit_is_in_fp:
          continue
        if evalue < best_sp_evalue:
          best_sp_hit = sp_hit_id
          best_sp_evalue = evalue
    if best_sp_hit != None:
      fp_str = '\t'
      if best_sp_hit_is_in_fp:
        fp_str = '*\t' + ';'.join(fp_hits[best_sp_hit])
      print "%s\t%s\t%s\t+\t%s" \
          % (best_sp_hit, best_sp_evalue, seed_id, fp_str)
    best_hits_for_organisms = {}
    best_evalues_for_organisms = {}
    best_hit_for_organism_is_in_fp = {}
    for organism in organisms:
      best_hits_for_organisms[organism] = ""
      best_evalues_for_organisms[organism] = 100
      best_hit_for_organism_is_in_fp[organism] = False
    for hit_line in hit_lines:
      hit_line = hit_line.rstrip()
      fields = hit_line.split()
      hit_id = fields[1].split('|')[2]
      is_fp_hit = hit_id in fp_hits and fp_hits[hit_id] != None
      organism = hit_id.split('_')[1]
      evalue = float(fields[10])
      if is_fp_hit:
        if not best_hit_for_organism_is_in_fp[organism] \
            or evalue < best_evalues_for_organisms[organism]:
          best_hit_for_organism_is_in_fp[organism] = True
          best_hits_for_organisms[organism] = hit_id
          best_evalues_for_organisms[organism] = evalue
      else:
        if evalue > evalue_cutoff or best_hit_for_organism_is_in_fp[organism]:
          continue
        if evalue < best_evalues_for_organisms[organism]:
          best_hits_for_organisms[organism] = hit_id
          best_evalues_for_organisms[organism] = evalue
    for organism in organisms:
      best_hit_id = best_hits_for_organisms[organism]
      if best_hit_id != "" and best_hit_id != best_sp_hit:
        best_evalue = best_evalues_for_organisms[organism]
        sp_str = ""
        if best_hit_id in sp_hit_set:
          sp_str = '+'
        fp_str = '\t'
        if best_hit_for_organism_is_in_fp[organism]:
          fp_str = '*\t' + ';'.join(fp_hits[best_hit_id])
        print "%s\t%s\t%s\t%s\t%s" \
            % (best_hit_id, best_evalue, seed_id, sp_str, fp_str)
    os.chdir(cur_dir)
  sp_pickle_f = open("SwissProtHits.pkl", "w")
  cPickle.dump(sp_hit_set, sp_pickle_f)
  sp_pickle_f.close()
  print "%d total SwissProt hits" % len(sp_hit_set)

if __name__ == '__main__':
  main()
