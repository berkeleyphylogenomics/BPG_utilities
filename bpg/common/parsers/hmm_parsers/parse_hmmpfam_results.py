#!/usr/bin/env python

import os
import re
from optparse import OptionParser


MAX_EVALUE = 0.001
MIN_Q_MATCH_LENGTH = 1
MIN_HMM_MATCH_LENGTH = 1


class DomainMatch:

  def __init__(self, name, evalue, qstart, qend, hstart, hend):
    self.name = name
    self.evalue = float(evalue)
    self.qstart = int(qstart)
    self.qend = int(qend)
    self.qlength = self.qend + 1 - self.qstart
    self.hstart = int(hstart)
    self.hend = int(hend)
    self.hlength = self.hend + 1 - self.hstart

  def __unicode__(self):
    return u"%s: %s [%d, %d]: %d" % (self.name, self.evalue, self.qstart,

        self.qend, self.qlength)
  def __str__(self):
    return self.__unicode__()


class ParsedResults:
  def __init__(self, query_sequence, domains):
    self.query_sequence = query_sequence
    self.domains = domains


def parse(input_file, max_evalue,
          min_query_match_length,
          min_hmm_match_length):

  f = open(input_file)
  lines = [line.rstrip() for line in f.readlines()]
  f.close()

  query_sequence_re = re.compile("Query sequence: (.*)")
  parsed_for_domains_re = re.compile("Parsed for domains:")

  SENTINEL = r'//'
  SEEKING_QUERY, SEEKING_PARSED, SKIPPING_MODEL_LINE, SKIPPING_RULER, \
    PARSING_DOMAINS = xrange(5)
  state = SEEKING_QUERY

  def advance_state(state):
    return (state + 1) % 5

  def cmp_domain_matches(match0, match1):
    result = cmp(match0.qstart, match1.qstart)
    if result == 0:
      result = cmp(match0.evalue, match1.evalue)
    if result == 0:
      result = cmp(match1.qlength, match0.qlength)
    if result == 0:
      result = cmp(match1.hlength, match0.hlength)
    if result == 0:
      result = cmp(match0.name, match1.name)
    if result == 0: # should never get here
      result = cmp(match0.hstart, match1.hstart)
    return result

  query_sequence = ""
  domains = []

  for line in lines:
    if state == SEEKING_QUERY:
      m = query_sequence_re.match(line)
      if m:
        query_sequence = m.group(1)
        state = advance_state(state)

    elif state == SEEKING_PARSED:
      if parsed_for_domains_re.match(line):
        state = advance_state(state)

    elif state == SKIPPING_MODEL_LINE or state == SKIPPING_RULER:
      state = advance_state(state)

    elif state == PARSING_DOMAINS:
      if len(line) == 0 or line == SENTINEL:
        domains.sort(cmp_domain_matches)

        return ParsedResults(query_sequence, domains)

        query_sequence = ""
        domains = []
        state = advance_state(state)
      else:
        if line.strip() == "[no hits above thresholds]":
          query_sequence = ""
          domains = []
          state = advance_state(state)
        else:
          name, count, qstart, qend, dots, hstart, hend, brackets, score,\
              evalue = line.split()
          domain_match = DomainMatch(name, evalue, qstart, qend, hstart, hend)

          if domain_match.evalue <= max_evalue and \
              domain_match.qlength >= min_query_match_length and \
              domain_match.hlength >= min_hmm_match_length:
            domains.append(domain_match)


def main():
  # parse command line options
  usage = "%prog [options] hmmpfam_results_file"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-e", "--max_evalue", dest="max_evalue",
        default=MAX_EVALUE, help="maximum e-value to include")
  opt_parser.add_option("-q", "--min_query_match_length",
        dest="min_query_match_length",
        default=MIN_Q_MATCH_LENGTH,
        help="minimum # of matching residues in the query")
  opt_parser.add_option("-d", "--min_hmm_match_length",
        dest="min_hmm_match_length",
        default=MIN_HMM_MATCH_LENGTH,
        help="minimum # of matching nodes in the hmm")

  (options, args) = opt_parser.parse_args()

  # verify input for errors
  if len(args) != 1:
    opt_parser.error("Incorrect number of arguments")
  input_file = args[0]
  if not os.path.exists(input_file):
    opt_parser.error("HMMPFAM results file %s does not exist" % input_file)

  # Validate max_evalue
  try:
    max_evalue = float(options.max_evalue)
  except ValueError:
    opt_parser.error("--max_evalue must be a nonnegative real")
  if max_evalue < 0:
    opt_parser.error("--max_evalue must be a nonnegative real")

  # Validate min_query_match_length
  try:
    min_query_match_length = int(options.min_query_match_length)
  except ValueError:
    opt_parser.error("--min_query_match_length must be a positive number")

  if min_query_match_length <= 0:
    opt_parser.error("--min_query_match_length must be a positive number")

  # Validate min_hmm_match_length
  try:
    min_hmm_match_length = int(options.min_hmm_match_length)
  except ValueError:
    opt_parser.error("--min_hmm_match_length must be a positive number")

  if min_hmm_match_length <= 0:
    opt_parser.error("--min_hmm_match_length must be a positive number")

  result = parse(input_file, max_evalue, min_query_match_length,
                  min_hmm_match_length)

  # Print results
  print "Query,Domain,E-value,QueryMatchStart,QueryMatchEnd," \
        + "HMMMatchStart,HMMMatchEnd"
  for domain in result.domains:
    print "%s,%s,%g,%d,%d,%d,%d" % (result.query_sequence, domain.name,
                                    domain.evalue,
                                    domain.qstart, domain.qend,
                                    domain.hstart, domain.hend)

if __name__ == '__main__':
  main()
