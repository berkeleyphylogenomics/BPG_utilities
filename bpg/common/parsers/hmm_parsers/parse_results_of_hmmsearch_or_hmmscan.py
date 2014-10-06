#!/usr/bin/env python
"""Parser to interpret the HMMER 3.0b2 output file

This parser will, in essence, return a python object representation of the
data given by HMMER 3.0. Currently, this file is written for HMMER 3.0b2
format. And, if HMMER changes its output format when it moves from Beta to
production, this parser may not function properly and will need to be
updated
"""

import os, sys, string, re
from optparse import OptionParser

MAX_EVALUE = 0.001
MIN_Q_MATCH_LENGTH = 1
MIN_HMM_MATCH_LENGTH = 1


NUM_STATES = 13
HEADER, SEEKING_COMPLETE_SEQUENCE_SCORES, \
     SEEKING_SEPARATOR, SCORE_INCLUSION, SCORE_EXCLUSION, \
     SEEKING_HIT_HEADER, \
     SEEKING_HEADER_LINE, SEEKING_HEADER_SEPARATOR, SEEKING_MATCH_INFO, \
     SEEKING_ALIGNMENTS_FOR_EACH_DOMAIN, SEEKING_DOMAIN_MATCH, \
     SEEKING_ALIGNED_SEQUENCE, INTERNAL_PIPELINE \
     = xrange(NUM_STATES)
name_of_state = {}
name_of_state[HEADER] = 'HEADER'
name_of_state[SEEKING_COMPLETE_SEQUENCE_SCORES] = 'SEEKING_COMPLETE_SEQUENCE_SCORES'
name_of_state[SEEKING_SEPARATOR] = 'SEEKING_SEPARATOR'
name_of_state[SCORE_INCLUSION] = 'SCORE_INCLUSION'
name_of_state[SCORE_EXCLUSION] = 'SCORE_EXCLUSION'
name_of_state[SEEKING_HIT_HEADER] = 'SEEKING_HIT_HEADER'
name_of_state[SEEKING_HEADER_LINE] = 'SEEKING_HEADER_LINE'
name_of_state[SEEKING_HEADER_SEPARATOR] = 'SEEKING_HEADER_SEPARATOR'
name_of_state[SEEKING_MATCH_INFO] = 'SEEKING_MATCH_INFO'
name_of_state[SEEKING_ALIGNMENTS_FOR_EACH_DOMAIN] = 'SEEKING_ALIGNMENTS_FOR_EACH_DOMAIN'
name_of_state[SEEKING_DOMAIN_MATCH] = 'SEEKING_DOMAIN_MATCH'
name_of_state[SEEKING_ALIGNED_SEQUENCE] = 'SEEKING_ALIGNED_SEQUENCE'
name_of_state[INTERNAL_PIPELINE] = 'INTERNAL_PIPELINE'

NUM_ALIGNMENT_STATES = 4
SEEKING_DOMAIN_ALIGNMENT_LINE, SEEKING_CONSENSUS_LINE, \
    SEEKING_SEQUENCE_ALIGNMENT_LINE, SEEKING_PP_LINE \
    = xrange(NUM_ALIGNMENT_STATES)

query_re = re.compile("Query: +(\S+) +")
complete_sequence_scores_re = re.compile("Scores for complete sequence")
inclusion_threshold_re = re.compile("-{6} inclusion threshold -{6}")
domain_and_alignment_re \
    = re.compile("Domain( and alignment)? annotation for each")
full_sequence_re = re.compile("full sequence")
evalue_re = re.compile("E-value")
cevalue_re = re.compile("c-Evalue")
separator_re = re.compile("-------")

hit_header_re = re.compile("^>>")
alignments_for_each_domain_re = re.compile("Alignments for each domain")
domain_match_re = re.compile("== domain ([0-9]*)")
internal_pipeline_re = re.compile("Internal pipeline statistics summary:")
sentinel_re = re.compile(r'//')
number_re = re.compile('\d+')


class HMMSearchOrScanResults:
    """Container of HitResult objects.

    This is basically a dictionary of HitResult objects.  It will
    only add matches that are significant.
    """
    def __init__(self, max_evalue, min_sequence_match_length,
                  min_hmm_match_length):
      self.max_evalue = max_evalue
      self.min_hmm_match_length = min_sequence_match_length
      self.min_sequence_match_length = min_sequence_match_length
      self.hit_result_of_name_of_query = {}

    def start_query(self, query):
      self.hit_result_of_name_of_query[query] = {}

    def add(self, query, name, hit_result):
      self.hit_result_of_name_of_query[query][name] = hit_result

    def add_match(self, query, name, number_str, inclusion_str, bit_score_str, 
                    i_evalue_str, hmm_from_str, hmm_to_str, 
                    seq_from_str, seq_to_str, match_type):
      if query not in self.hit_result_of_name_of_query:
        return False
      if name not in self.hit_result_of_name_of_query[query]:
        return False
      number = int(number_str)
      inclusion = (inclusion_str == '!')
      bit_score = float(bit_score_str)
      i_evalue = float(i_evalue_str)
      hmm_from = int(hmm_from_str)
      hmm_to = int(hmm_to_str)
      seq_from = int(seq_from_str)
      seq_to = int(seq_to_str)
      if i_evalue <= self.max_evalue and \
          (hmm_to - hmm_from >= self.min_hmm_match_length) and \
          (seq_to - seq_from >= self.min_sequence_match_length):
        self.hit_result_of_name_of_query[query][name].add_match(number, 
                      inclusion,
                      bit_score, i_evalue, hmm_from, hmm_to, seq_from, seq_to,
                      match_type)
        return True
      return False

    def update_match_aligned_chars(self, query, name, number, aligned_hit,
                                    num_aligned_chars):
      if query in self.hit_result_of_name_of_query and \
          name in self.hit_result_of_name_of_query[query]:
        self.hit_result_of_name_of_query[query][
            name].update_match_aligned_chars(
                                      number, aligned_hit, num_aligned_chars)
      elif query not in self.hit_result_of_name_of_query:
        sys.stderr.write("query %s missing\n" % query)

    def remove(self, query, name):
      if query in self.hit_result_of_name_of_query and \
          name in self.hit_result_of_name_of_query[query]:
        del self.hit_result_of_name_of_query[query][name]

    def __unicode__(self):
      return '; '.join(['Query %s: %s' % (query,
                        ''.join([result.__str__() for result in
                        self.hit_result_of_name_of_query[query].values()]))
                        for query in self.hit_result_of_name_of_query.keys()])

    def __str__(self):
      return self.__unicode__()

class MatchResult:
  """Single match of a domain to a sequence."""
  def __init__(self, inclusion, bit_score, i_evalue, hmm_from, hmm_to,
              seq_from, seq_to, match_type):
      self.inclusion = inclusion
      self.bit_score = bit_score
      self.i_evalue = i_evalue
      self.hmm_from = hmm_from
      self.hmm_to = hmm_to
      self.seq_from = seq_from
      self.seq_to = seq_to
      self.match_type = match_type
      self.num_aligned_chars = None
      self.aligned_hit = None

  def set_aligned_hit(self, aligned_hit):
      self.aligned_hit = aligned_hit

  def set_num_aligned_chars(self, num_aligned_chars):
      self.num_aligned_chars = num_aligned_chars

  def __unicode__(self):
    if self.inclusion:
      inclusion_str = 'Y'
    else:
      inclusion_str = 'N'
    result =  "Include: %s " % str(self.inclusion) \
            + "Bit score: %g " % self.bit_score \
            + "E-value: %g " % self.i_evalue \
            + "HMM Start: %d " % self.hmm_from \
            + "HMM End: %d " % self.hmm_to \
            + "Sequence Start: %d " % self.seq_from \
            + "Sequence End: %d " % self.seq_to \
            + "# Aligned Chars: %d " % self.num_aligned_chars \
            + "Match Type: %s" % self.match_type
    return result

  def __str__(self):
    return self.__unicode__()

class HitResult:
    """Single object from hmmsearch output file summary section."""

    def __init__(self, full_evalue, full_score, full_bias,
                       best_evalue, best_score, best_bias,
                       exp, n, name, description,
                       inclusion=True):

        self.full_evalue = float(full_evalue)
        self.full_score = float(full_score)
        self.full_bias = float(full_bias)

        self.best_evalue = float(best_evalue)
        self.best_score = float(best_score)
        self.best_bias = float(best_bias)

        self.exp = float(exp)
        self.n = int(n)
        self.name = name
        self.description = description
        self.inclusion = inclusion
        self.matches = {}

    def add_match(self, number, inclusion, bit_score, i_evalue,
                  hmm_from, hmm_to, seq_from, seq_to, match_type):
      match_result = MatchResult(inclusion, bit_score, i_evalue, 
                                hmm_from, hmm_to, seq_from, seq_to, match_type)
      self.matches[number] = match_result

    def update_match_aligned_chars(self, number, aligned_hit,
                                    num_aligned_chars):
      if number in self.matches:
        self.matches[number].set_num_aligned_chars(num_aligned_chars)
        self.matches[number].set_aligned_hit(aligned_hit)
        return True
      else:
        return False

    def __unicode__(self):
        return u"%s: %s\n" % (self.name, self.description) \
              + ''.join([u"%d) %s\n" % (n, self.matches[n].__unicode__())
                            for n in self.matches.keys()])

    def __str__(self):
        return self.__unicode__()

def process_score_line(line, inclusion=False):

    fields = line.split()
    if len(fields) >= 9:
        full_evalue, full_score, full_bias, best_evalue, best_score, \
        best_bias, exp, n, name = fields[:9]
        description = " ".join(fields[9:])

        result = HitResult(full_evalue, full_score, full_bias,
                                best_evalue, best_score, best_bias,
                                exp, n, name, description, True)

        return name, result 
    return (None, None)


def parse(input_file, max_evalue, min_sequence_match_length, 
          min_hmm_match_length):
    hmm_search_or_scan_results = HMMSearchOrScanResults(max_evalue, 
                                          min_sequence_match_length,
                                          min_hmm_match_length)

    f = open(input_file)
    lines = [line.strip() for line in f.readlines()]
    f.close()

    trivial_translation = string.maketrans('', '')
    dotdashlowercase = '.-' + string.lowercase

    state = HEADER
    alignment_state = SEEKING_DOMAIN_ALIGNMENT_LINE


    query = ''
    name = ''
    aligned_chars = ''
    domain_number = None
    num_queries = 0
    def advance_state(state):
        new_state = (state + 1) % NUM_STATES
        return new_state
    def advance_alignment_state(alignment_state):
        new_alignment_state = (alignment_state + 1) % NUM_ALIGNMENT_STATES
        return new_alignment_state

    for line in lines:

        if state == HEADER:
            m = query_re.search(line)
            if m:
                query = m.group(1)
                hmm_search_or_scan_results.start_query(query)
                num_queries += 1
                state = advance_state(state)

        if state == SEEKING_COMPLETE_SEQUENCE_SCORES:
            if complete_sequence_scores_re.search(line):
                state = advance_state(state)

        elif state == SEEKING_SEPARATOR:
            if separator_re.search(line):
                state = advance_state(state)

        elif state == SCORE_INCLUSION or state == SCORE_EXCLUSION:
            if domain_and_alignment_re.search(line): 
                state = SEEKING_HIT_HEADER
            elif inclusion_threshold_re.search(line):
                state = advance_state(state)
            elif not (full_sequence_re.search(line) or \
                     evalue_re.search(line) or \
                     separator_re.search(line)):
              name, hit_result = process_score_line(line,
                                            inclusion=(state==SCORE_INCLUSION))
              if hit_result:
                  if hit_result.best_evalue <= max_evalue:
                    hmm_search_or_scan_results.add(query, name, hit_result)

        elif state == SEEKING_HIT_HEADER:
            if internal_pipeline_re.search(line):
                state = INTERNAL_PIPELINE
            elif hit_header_re.match(line):
                name = line[2:].strip().split()[0]
                state = advance_state(state)

        elif state == SEEKING_HEADER_LINE:
            if cevalue_re.search(line):
              state = advance_state(state)
            elif internal_pipeline_re.match(line):
              state = INTERNAL_PIPELINE
            elif hit_header_re.match(line):
              name = line[2:].strip().split()[0]
              # The state is still SEEKING_HEADER_LINE
              

        elif state == SEEKING_HEADER_SEPARATOR:
            if separator_re.search(line):
              state = advance_state(state)

# Example:
# >> 3hhm_B  mol:protein length:373  niSH2 p85alpha
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#    1 !   43.9   0.0   4.3e-13   1.9e-10       1      86 []      13      87 ..      13      87 .. 0.91

        elif state == SEEKING_MATCH_INFO:
            if len(line) == 0:
              state = advance_state(state)
            else:
              fields = line.strip().split()
              hmm_search_or_scan_results.add_match(query, name, fields[0], 
                                  fields[1], fields[2], fields[5], 
                                  fields[6], fields[7], fields[9], 
                                  fields[10], fields[8])
              
        elif state == SEEKING_ALIGNMENTS_FOR_EACH_DOMAIN:
            if alignments_for_each_domain_re.search(line):
              state = advance_state(state)

        elif state == SEEKING_DOMAIN_MATCH:
            m = domain_match_re.search(line)
            if m:
              domain_number = int(m.group(1))
              aligned_sequence = ''
              alignment_state = SEEKING_DOMAIN_ALIGNMENT_LINE
              state = advance_state(state)

        elif state == SEEKING_ALIGNED_SEQUENCE:
            if internal_pipeline_re.search(line):
              if domain_number:
                aligned_chars = aligned_sequence.translate(
                                                         trivial_translation, 
                                                          dotdashlowercase)
                hmm_search_or_scan_results.update_match_aligned_chars(query, 
                                                            name, 
                                                            domain_number,
                                                            aligned_sequence,
                                                            len(aligned_chars))
              aligned_sequence = ''
              domain_number = None
              state = advance_state(state)
              continue
            if len(line) > 0:
              # Keep looping through the 5-line alignment blocks
              fields = line.split()
              if alignment_state == SEEKING_CONSENSUS_LINE:
                alignment_state = advance_alignment_state(alignment_state)
              elif alignment_state == SEEKING_SEQUENCE_ALIGNMENT_LINE:
                aligned_sequence += fields[2]
                alignment_state = advance_alignment_state(alignment_state)
              elif alignment_state == SEEKING_PP_LINE:
                alignment_state = advance_alignment_state(alignment_state)
              elif alignment_state == SEEKING_DOMAIN_ALIGNMENT_LINE:
                if fields[1] == 'CS':
                  continue
                elif len(fields) == 4 and number_re.match(fields[1]) and \
                    number_re.match(fields[3]):
                  alignment_state = advance_alignment_state(alignment_state)
                elif hit_header_re.match(line):
                  aligned_chars = aligned_sequence.translate(
                                                           trivial_translation, 
                                                            dotdashlowercase)
                  hmm_search_or_scan_results.update_match_aligned_chars(query,
                                                          name,
                                                          domain_number,
                                                          aligned_sequence,
                                                          len(aligned_chars))
                  name = line[2:].strip().split()[0]
                  aligned_sequence = ''
                  state = SEEKING_HEADER_LINE
                  # The alignment_state is still SEEKING_DOMAIN_ALIGNMENT_LINE
                else:
                  m = domain_match_re.search(line)
                  if m:
                    aligned_chars = aligned_sequence.translate(
                                                          trivial_translation, 
                                                          dotdashlowercase)
                    hmm_search_or_scan_results.update_match_aligned_chars(
                                                            query,
                                                            name,
                                                            domain_number,
                                                            aligned_sequence,
                                                            len(aligned_chars))
                    domain_number = int(m.group(1))
                    aligned_sequence = ''
                    # The state is again SEEKING_ALIGNED_SEQUENCE
                    # The alignment_state is still 
                    # SEEKING_DOMAIN_ALIGNMENT_LINE
            elif alignment_state == SEEKING_CONSENSUS_LINE:
              # The consensus line may be blank
              alignment_state = advance_alignment_state(alignment_state)
          
        elif state == INTERNAL_PIPELINE:
            if sentinel_re.search(line):
              # Loop around to the beginning for the next query
              state = advance_state(state)

    # When the parser is used against files where No hits are detected (as below)
    # then the domain_number will not be defined and the update_match_aligned_chars
    # fails. This was fixed by setting domain_number to None by default and only
    # calling this next function if domain_number is not None.
    #
    # grj -- 2010 Apr 5
    #
    # Scores for complete sequences (score includes all domains):
    #   --- full sequence ---   --- best 1 domain ---    -#dom-
    #    E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
    #    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
    #
    #   [No hits detected that satisfy reporting thresholds]
    #
    #
    # Domain and alignment annotation for each sequence:
    #
    #   [No hits detected that satisfy reporting thresholds]

    if domain_number:
        aligned_chars = aligned_sequence.translate(trivial_translation, 
                                                  dotdashlowercase)
        hmm_search_or_scan_results.update_match_aligned_chars(query, name, 
                                                      domain_number,
                                                      aligned_sequence,
                                                      len(aligned_chars))
    return hmm_search_or_scan_results

def main():
    # parse command line options
    usage = "%prog <hmmpfam_results_file>"
    opt_parser = OptionParser(usage=usage)
    opt_parser.add_option("-e", "--max_evalue", dest="max_evalue",
          default=MAX_EVALUE, help="maximum e-value to include")
    opt_parser.add_option("-q", "--min_sequence_match_length",
          dest="min_sequence_match_length",
          default=MIN_Q_MATCH_LENGTH,
          help="minimum # of matching residues in the sequence")
    opt_parser.add_option("-d", "--min_hmm_match_length",
          dest="min_hmm_match_length",
          default=MIN_HMM_MATCH_LENGTH,
          help="minimum # of matching nodes in the hmm")

    (options, args) = opt_parser.parse_args()

    # verify input for errors
    if len(args) != 1:
        opt_parser.error("Include filename to parse")

    input_file = args[0]

    # Validate max_evalue
    try:
      max_evalue = float(options.max_evalue)
    except ValueError:
      opt_parser.error("--max_evalue must be a nonnegative real")
    if max_evalue < 0:
      opt_parser.error("--max_evalue must be a nonnegative real")

    # Validate min_sequence_match_length
    try:
      min_sequence_match_length = int(options.min_sequence_match_length)
    except ValueError:
      opt_parser.error("--min_sequence_match_length must be a positive number")

    if min_sequence_match_length <= 0:
      opt_parser.error("--min_sequence_match_length must be a positive number")

    # Validate min_hmm_match_length
    try:
      min_hmm_match_length = int(options.min_hmm_match_length)
    except ValueError:
      opt_parser.error("--min_hmm_match_length must be a positive number")

    if min_hmm_match_length <= 0:
      opt_parser.error("--min_hmm_match_length must be a positive number")

    # Parse file and print results
    hmm_search_or_scan_results = parse(input_file, max_evalue, 
                              min_sequence_match_length, min_hmm_match_length)
    sys.stdout.write("QueryId,HitId,MatchNumber,Include?,E-value,BitScore,")
    sys.stdout.write("HMMStart,HMMEnd,SequenceStart,SequenceEnd,")
    sys.stdout.write("NumAlignedChars,MatchType\n")
    for query in hmm_search_or_scan_results.hit_result_of_name_of_query:
      for name in hmm_search_or_scan_results.hit_result_of_name_of_query[query]:
        for number in hmm_search_or_scan_results.hit_result_of_name_of_query[
                                                          query][name].matches:
          match_result = hmm_search_or_scan_results.hit_result_of_name_of_query[
                                                  query][name].matches[number]
          sys.stdout.write("%s,%s,%d," % (query, name, number))
          if match_result.inclusion:
            sys.stdout.write("Y,")
          else:
            sys.stdout.write("N,")
          sys.stdout.write("%g,%g,%d,%d,%d,%d,%d,%s\n"
            % (match_result.i_evalue, match_result.bit_score,
              match_result.hmm_from, match_result.hmm_to,
              match_result.seq_from, match_result.seq_to,
              match_result.num_aligned_chars, match_result.match_type))
          sys.stdout.flush()

if __name__ == '__main__':
    main()
