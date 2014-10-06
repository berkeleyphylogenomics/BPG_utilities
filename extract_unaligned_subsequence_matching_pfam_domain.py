#!/usr/bin/env python

from bpg.common.parsers.hmm_parsers \
    import parse_results_of_hmmsearch_or_hmmscan
from optparse import OptionParser
import string

def main():
  usage = "%prog [options] results_of_hmmscan_of_single_seq_vs_pfam"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')

  hmmscan_results = parse_results_of_hmmsearch_or_hmmscan.parse(args[0], 0.001,
                                                            1, 1)
  uppercase_translation = string.maketrans(string.lowercase, string.uppercase)
  dotdash = '.-'
  #:print "Query,Pfam,Match#,SeqFrom,SeqTo,UnalignedHit"
  for query in hmmscan_results.hit_result_of_name_of_query:
    for pfam_name in hmmscan_results.hit_result_of_name_of_query[query]:
      hit_result = hmmscan_results.hit_result_of_name_of_query[query][pfam_name]
      for match_number in hit_result.matches:
        match_result = hit_result.matches[match_number]
        unaligned_hit = match_result.aligned_hit.translate(
                                uppercase_translation, dotdash)
        print ">%s,%s,%s,%s,%s\n%s" % (query,pfam_name, match_number, 
                                  match_result.seq_from, match_result.seq_to, 
                                  unaligned_hit)

if __name__ == '__main__':
  main()
