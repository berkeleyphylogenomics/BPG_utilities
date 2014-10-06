#!/usr/bin/env python

from optparse import OptionParser
import os
import re
import string
import fcntl

from Bio import Seq, SeqIO

from bpg.common import BPGPWID, fasta
from bpg.common.get_blastable_database_parameters import get_blastable_database_parameters

def minimum_coverage(sequence_length):
  if sequence_length < 100:     
    return 0.60 + 0.0
  elif sequence_length < 200:
    return 0.60 + 0.05
  elif sequence_length < 250:
    return 0.60 + 0.10
  elif sequence_length < 300:
    return 0.60 + 0.13
  elif sequence_length < 350:
    return 0.60 + 0.15
  elif sequence_length < 400:
    return 0.60 + 0.18
  elif sequence_length < 450:
    return 0.60 + 0.20
  elif sequence_length < 500:
    return 0.60 + 0.23
  else:           # >= 500
    return 0.60 + 0.25 


trivial_translation = string.maketrans('', '')
contiguous_gap = re.compile('-+')
contiguous_lowercase = re.compile('[a-z]+')
dotdash = '.-'
dotlowercase = '.' + string.lowercase
dotdashlowercase = dotdash + string.lowercase
def test_coverage(hit_seq, hit_length, seed_length):
  upper_only = hit_seq.translate(trivial_translation, dotdashlowercase)
  num_aligned = len(upper_only)
  # Don't check the proportion of the hit covered, since this is a glocal
  # search.
  if num_aligned < minimum_coverage(seed_length) * seed_length:
    return False
  cols_only = hit_seq.translate(trivial_translation, dotlowercase)
  gappy_regions = contiguous_gap.findall(cols_only)
  if len(gappy_regions) > 0:
    longest_gap = max(len(r) for r in gappy_regions)
    if longest_gap >= 40:
      return False
  letters_only = hit_seq.translate(trivial_translation, dotdash)
  for lowercase_region_match in contiguous_lowercase.finditer(letters_only):
    # Don't check the initial or final lowercase region, since this is a glocal
    # search.
    if lowercase_region_match.start > 0 and \
        lowercase_region_match.end < len(letters_only):
      if lowercase_region_match.end - lowercase_region_match.start >= 40:
        return False
  return True

def main():
  # parse command line options
  usage = "%prog [options] pfam_subdir_to_check"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-n", "--num_hmm_iterations", 
              dest="num_hmm_iterations",
              default=3, help="Number of HMM iterations")
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  try:
    num_hmm_iterations = int(options.num_hmm_iterations)
  except ValueError:
    opt_parser.error("--num_hmm_iterations must be a number")
  pfam_seed_to_check = args[0]
  os.chdir('/clusterfs/ohana/bpg/pfam')
  pfam_subdir, seed_file = os.path.split(pfam_seed_to_check)
  if pfam_subdir.find('subfam') >= 0:
    universe_subdir, seed_subdir = os.path.split(pfam_subdir)
    relative_universe_subdir = '..'
  else:
    universe_subdir, seed_subdir = pfam_subdir, '.'
    relative_universe_subdir = '.'
  os.chdir(universe_subdir)
  universe_fa = 'universe-01.fa'
  if not os.path.exists('%s.pin' % universe_fa):
    cmd = 'formatdb -o F -i %s -l %s.log >& %s.err' \
      % (universe_fa, universe_fa, universe_fa)
    print cmd
    os.system(cmd)
  os.chdir(seed_subdir)
  universe_fa = os.path.join(relative_universe_subdir, universe_fa)
  seed_handle = open(seed_file, "rU")
  seed_key, seed_record = SeqIO.to_dict(SeqIO.parse(seed_handle, 
                                                          "fasta")).popitem()
  seed_handle.close()
  seed_id = seed_key
  # Replace / in the header by %2f
  seed_id = seed_id.replace('/','_')
  seed_length = len(seed_record.seq)
  universe_handle = open(universe_fa, "rU")
  universe = SeqIO.to_dict(SeqIO.parse(universe_handle, "fasta"))
  universe_handle.close()
  key_of_id = {}
  for key in universe.keys():
    key_of_id[key] = key
  max_blast_evalue = 1e6
  evalue_cutoff = 1e-5
  if (seed_length <= 100):
    evalue_cutoff = 10 * evalue_cutoff
  if (seed_length < 65):
    evalue_cutoff = 10 * evalue_cutoff
  if (seed_length < 15):
    evalue_cutoff = 10 * evalue_cutoff
  blast_result_filename = "%s_qa_big_universe.blast_results" \
                        % seed_key.split('/')[0]
  blast_exe = "blastall -p blastp"
  num_seqs, uniprot_db_length, long_len = get_blastable_database_parameters()
  blast_exe = blast_exe + " -z %d " % uniprot_db_length
  cmd = blast_exe \
    + " -i \"%s\" -o \"%s\" -m 9 -I T -b 4000 -v 4000 -e %f -d \"%s\"" \
    % (seed_file, blast_result_filename, max_blast_evalue, universe_fa)
  print cmd
  os.system(cmd)
  f = open(blast_result_filename, "rU")
  lines = f.readlines()
  f.close()
  hits = {}
  loose_hits = {}
  evalues = {}
  for line in lines:
    line = line.rstrip()
    if line[0] == '#':
      continue
    fields = line.split()
    hit_id = fields[1]
    if hit_id in loose_hits.keys():
      continue
    hit_key = key_of_id[hit_id]
    hit_record = universe[hit_key]
    hit_length = len(hit_record.seq)
    # Escape commas - align2model converts them into spaces
    hit_record.id = hit_record.id.replace(',','%2c')
    if hit_length >= minimum_coverage(seed_length) * seed_length:
      loose_hits[hit_id] = hit_record
      try:
        e_value = float(fields[10])
        if e_value <= evalue_cutoff:
            hits[hit_id] = hit_record
            evalues[hit_id] = e_value
      except ValueError:
        pass
  blast_homologs_filename = "%s_qa_blast_homologs.fa" % seed_id
  f = open(blast_homologs_filename, "w")
  SeqIO.write(hits.values(), f, "fasta")
  f.close()
  loose_blast_homologs_filename = "%s_qa_loose_blast_homologs.fa" % seed_id
  f = open(loose_blast_homologs_filename, "w")
  SeqIO.write(loose_hits.values(), f, "fasta")
  f.close()
  cmd = "w0.5 \"%s.fa\" \"%s.mod\" >& qa_w0.5.out" % (seed_id, seed_id)
  while not os.path.exists("%s.mod" % seed_id):
    print cmd
    os.system(cmd)
  cmd = "align2model \"%s_qa_iter0_homologs\" -i \"%s.mod\" -db \"%s\" " \
        % (seed_id, seed_id, blast_homologs_filename) \
        + " -sw 2 -adpstyle 5 >& qa_align2model.out"
  print cmd
  os.system(cmd)
  
  if len(hits) > 0:
    prev_iter = 0
    cur_iter = 1
    while cur_iter <= num_hmm_iterations:
      hmm_file = "%s_qa_iter%d_homologs.mod" % (seed_id, prev_iter)
      cmd = "w0.5 \"%s_qa_iter%d_homologs.a2m\" \"%s\" " \
            % (seed_id, prev_iter, hmm_file) \
            + " >& w0.5.qa_iter%d.out" % prev_iter
      print cmd
      while not os.path.exists(hmm_file):
        os.system(cmd)
      cmd = "hmmscore %s_qa_iter%d -db \"%s\" -i \"%s\" " \
            % (seed_id, cur_iter, loose_blast_homologs_filename, hmm_file) 
      cmd = cmd + "-dbsize 100000 -sw 2 > %s_qa_iter%d.hmmscore.out 2>&1" \
                  % (seed_id, cur_iter)
      os.system(cmd)
      distfile = "%s_qa_iter%d.dist" % (seed_id, cur_iter)
      f = open(distfile, "rU")
      dist_lines = f.readlines()
      f.close()
      cur_iter_hits = []
      for line in dist_lines:
        if line == "":
          break
        if line[0] != '%':
          # for every entry in the dist file,
          # get its sequence id (SAM only prints the id up to first space)

          s = string.split(line)
          seq_id = s[0]
          # Unescape commas
          seq_id = seq_id.replace('%2c',',')

          # get simple, reverse, and e-value scores

          e_value = string.atof(s[4])
          if e_value <= evalue_cutoff:
            hit_record = universe[seq_id]
            hit_length = len(hit_record.seq)
            if hit_length >= minimum_coverage(seed_length) * seed_length \
                and seed_length >= minimum_coverage(hit_length) * hit_length:
              print seq_id
              evalues[seq_id] = e_value
              cur_iter_hits.append(hit_record)
      cur_iter_homologs_filename = "%s_qa_iter%d_homologs.fa" \
                                    % (seed_id, cur_iter)
      f = open(cur_iter_homologs_filename, "w")
      SeqIO.write(cur_iter_hits, f, "fasta")
      f.close()
      cmd = "align2model %s_qa_iter%d_homologs -i \"%s\" -db \"%s\" " \
            % (seed_id, cur_iter, hmm_file, cur_iter_homologs_filename) \
            + "-sw 2 -adpstyle 5 > align2model_qa_iter%d.out 2>&1" % cur_iter
      os.system(cmd)
      prev_iter = cur_iter
      cur_iter = cur_iter + 1
    cmd = "ln -s \"%s_qa_iter%d_homologs.a2m\" \"%s_qa_homologs.a2m\"" \
            % (seed_id, num_hmm_iterations, seed_id)
    print cmd
    os.system(cmd)
    cmd = "ln -s \"%s\" \"%s_qa_last.mod\"" % (hmm_file, seed_id)
    print cmd
    os.system(cmd)
  else:
    cmd = "ln -s \"%s_qa_iter0_homologs.a2m\" \"%s_qa_homologs.a2m\"" \
            % (seed_id, seed_id)
    print cmd
    os.system(cmd)
    cmd = "ln -s \"%s.mod\" \"%s_qa_last.mod\"" % (seed_id, seed_id)
    print cmd
    os.system(cmd)
  possible_homologs = fasta.ReadSequencesList("%s_qa_homologs.a2m" % seed_id)
  final_hits_with_pwids = [(1.0, seed_key)]
  seed_seq = seed_record.seq.tostring()
  print "Checking %d possible homologs" % len(possible_homologs)
  for hit_id, hit_seq in possible_homologs:
    # Unescape commas
    hit_id = hit_id.replace('%2c',',')
    hit_key = hit_id.rstrip().split(' ')[0]
    hit_length = len(universe[hit_key].seq)
    if test_coverage(hit_seq, hit_length, seed_length):
      if hit_key != seed_key:
        pwid = BPGPWID.pairwise_identity_KS_1(seed_seq, hit_seq)
        final_hits_with_pwids.append((pwid, hit_key))
  def compare_seqs(pwid_and_seq0, pwid_and_seq1):
    return cmp(pwid_and_seq0[0], pwid_and_seq1[0])
  final_hits_with_pwids.sort(compare_seqs, reverse=True)
  def record_of(key):
    if key == seed_key:
      return seed_record
    else:
      return universe[key]
  final_hits = [record_of(hit_key) for pwid, hit_key in final_hits_with_pwids]
  cluster_file_name = "%s_qa_cluster.fa" % seed_id
  f = open(cluster_file_name, "w")
  SeqIO.write(final_hits, f, "fasta")
  f.close()
  if os.path.exists('final_trimmed.fa'):
    family_handle = open('final_trimmed.fa', "rU")
    family = SeqIO.to_dict(SeqIO.parse(family_handle, "fasta"))
    family_handle.close()
  else:
    print "Family file final_trimmed.fa is missing."
    family = {}
  print "Family keys"
  for seq_id in family:
    print seq_id
  num_missing = 0
  num_included = 0
  f = open('qa_detailed.csv', 'w')
  f.write('SequenceId,PercentIdToSeed,EValueToHMM,PresentInFamily\n')
  for pwid, hit_key in final_hits_with_pwids:
    if hit_key == seed_key:
      f.write('%s,1.0,0.0,')
    else:
      f.write('%s,%0.3f,%g,' % (hit_key, pwid, evalues[hit_key]))
    if hit_key in family:
      f.write('Y\n')
      num_included += 1
    else:
      f.write('N\n')
      num_missing += 1
  f.close()
  f = open('/clusterfs/ohana/bpg/pfam/qa_summary.report', 'a')
  fcntl.lockf(f, fcntl.LOCK_EX)
  f.write('%s,%d,%d\n' % (pfam_subdir, num_included, num_missing))
  f.flush()
  fcntl.lockf(f, fcntl.LOCK_UN)
  f.close()

if __name__ == '__main__':
  main()
