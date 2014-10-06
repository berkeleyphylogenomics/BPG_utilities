#!/usr/bin/env python

import os, sys
from Bio import SeqIO
from bpg.common.parsers.hmm_parsers \
    import parse_results_of_hmmsearch_or_hmmscan

def main():
  dir = '/clusterfs/ohana/bpg/InitialMicrobialTargetSet/Pfams/'
  os.chdir(dir)
  f = open('initial_target_microbial_genomes.fa')
  target_seqs = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
  f.close()
  hmmscan_filename = "initial_target_microbial_genomes_vs_Pfam-A.hmmscan.out"
  hmmscan_results = parse_results_of_hmmsearch_or_hmmscan.parse(
                      hmmscan_filename, 0.001, 1, 1)
  for query in target_seqs.keys():
    query_residues = {}
    seq = target_seqs[query].seq.tostring()
    qlen = len(seq)
    for i in xrange(qlen):
      query_residues[i] = False
    for pfam_name in hmmscan_results.hit_result_of_name_of_query[query]:
      hit_result = hmmscan_results.hit_result_of_name_of_query[query][pfam_name]
      for match_number in hit_result.matches:
        match_result = hit_result.matches[match_number]
        i = match_result.seq_from - 1
        for j in xrange(len(match_result.aligned_hit)):
          if match_result.aligned_hit[j].isupper():
            query_residues[i] = True
          if match_result.aligned_hit[j].isalpha():
            i += 1
    range_end = -1
    range_start = -1
    for i in xrange(qlen):
      if query_residues[i]:
        if range_end >= 0:
          if range_end - range_start >= 30:
            sys.stdout.write('>%s/%s-%s\n' 
                              % (query, range_start+1, range_end+1))
            sys.stdout.write('%s\n' % (seq[range_start:range_end+1]))
          range_end = -1
          range_start = -1
      else:
        if range_start < 0:
          range_start = i
        range_end = i
    if range_end >= 0:
      if range_end - range_start >= 30:
        sys.stdout.write('>%s/%s-%s\n' % (query, range_start+1, range_end+1))
        sys.stdout.write('%s\n' % (seq[range_start:range_end+1]))

      

        

if __name__ == '__main__':
  main()
