#!/usr/bin/env python

from pfacts003.phog.orthologs import getOrthologQuerySet
from pfacts003.phog.views import PHOGRow
from pfacts003.phylofacts.models import OrthologTypes, preset_thresholds
from django.db import connection

def print_query_range(lower, upper):
  for query in connection.queries[lower:upper]:
    print "Time: %s" % query['time']
    print "Query: %s" % query['sql']

def main():
  threshold = preset_thresholds[OrthologTypes.PHOG_T_Tight]
  ortholog_type = OrthologTypes.PHOG_T_Tight
  (ownUniProt, phogs, best_phogs,
    orthologs, phog_of_ortholog, error) \
      = getOrthologQuerySet('OXYR_HUMAN', ortholog_type, threshold)
  print "Got %d tight orthologs" % len(orthologs)
  num_queries_for_tight_orthologs = len(connection.queries)
  cum_queries_for_tight_orthologs = len(connection.queries)
  print "%d queries were made to get tight orthologs" \
          % num_queries_for_tight_orthologs
  tight_orthologs_time = 0.0
  for query in connection.queries:
    tight_orthologs_time += float(query['time'])
  print "Total time for tight orthologs: %0.3f" % tight_orthologs_time
  phog_rows = [PHOGRow(phog, ortholog_type, threshold) for phog in best_phogs]
  print "Got %d tight PHOGRows" % len(phog_rows)
  num_queries_for_tight_phogrows \
    = len(connection.queries) - cum_queries_for_tight_orthologs
  cum_queries_for_tight_phogrows = len(connection.queries)
  print "%d queries were made to get tight PHOGRows" \
          % num_queries_for_tight_phogrows
  tight_phogrows_time = 0.0
  for query in connection.queries[cum_queries_for_tight_orthologs:]:
    tight_phogrows_time += float(query['time'])
  print "Total time for tight PHOGRows: %0.3f" % tight_phogrows_time
  for phog_row in phog_rows:
    print phog_row.accession, "%d sequences" \
      % phog_row.num_nonredundant_sequences
  threshold = 0.1
  ortholog_type = OrthologTypes.PHOG_T_Custom
  (ownUniProt, phogs, best_phogs,
    orthologs, phog_of_ortholog, error) \
      = getOrthologQuerySet('OXYR_HUMAN', ortholog_type, threshold)
  print "Got %d custom orthologs at threshold 0.1" % len(orthologs)
  num_queries_for_custom_orthologs \
    = len(connection.queries) - cum_queries_for_tight_phogrows
  cum_queries_for_custom_orthologs = len(connection.queries)
  print "%d queries were made to get custom orthologs" \
          % num_queries_for_custom_orthologs
  custom_orthologs_time = 0.0
  for query in connection.queries[cum_queries_for_tight_phogrows:]:
    custom_orthologs_time += float(query['time'])
  print "Total time for custom orthologs: %0.3f" % custom_orthologs_time
  phog_rows = [PHOGRow(phog, ortholog_type, threshold) for phog in best_phogs]
  print "Got %d custom PHOGRows" % len(phog_rows)
  num_queries_for_custom_phogrows \
    = len(connection.queries) - cum_queries_for_custom_orthologs
  cum_queries_for_custom_phogrows = len(connection.queries)
  print "%d queries were made to get custom PHOGRows" \
          % num_queries_for_custom_phogrows
  custom_phogrows_time = 0.0
  for query in connection.queries[cum_queries_for_custom_orthologs:]:
    custom_phogrows_time += float(query['time'])
  print "Total time for custom PHOGRows: %0.3f" % custom_phogrows_time
  for phog_row in phog_rows:
    print phog_row.accession, "%d sequences" \
      % phog_row.num_nonredundant_sequences
  threshold = preset_thresholds[OrthologTypes.PHOG_T_Loose]
  ortholog_type = OrthologTypes.PHOG_T_Loose
  (ownUniProt, phogs, best_phogs,
    orthologs, phog_of_ortholog, error) \
      = getOrthologQuerySet('OXYR_HUMAN', ortholog_type, threshold)
  print "Got %d loose orthologs" % len(orthologs)
  num_queries_for_loose_orthologs = \
      len(connection.queries) - cum_queries_for_custom_phogrows
  cum_queries_for_loose_orthologs = len(connection.queries)
  print "%d queries were made to get loose orthologs" \
          % num_queries_for_loose_orthologs
  loose_orthologs_time = 0.0
  for query in connection.queries[cum_queries_for_custom_phogrows:]:
    loose_orthologs_time += float(query['time'])
  print "Total time for loose orthologs: %0.3f" % loose_orthologs_time
  phog_rows = [PHOGRow(phog, ortholog_type, threshold) for phog in best_phogs]
  print "Got %d loose PHOGRows" % len(phog_rows)
  num_queries_for_loose_phogrows \
    = len(connection.queries) - cum_queries_for_loose_orthologs
  cum_queries_for_loose_phogrows = len(connection.queries)
  print "%d queries were made to get loose PHOGRows" \
          % num_queries_for_loose_phogrows
  loose_phogrows_time = 0.0
  for query in connection.queries[cum_queries_for_loose_orthologs:]:
    loose_phogrows_time += float(query['time'])
  print "Total time for loose PHOGRows: %0.3f" % loose_phogrows_time
  for phog_row in phog_rows:
    print phog_row.accession, "%d sequences" \
      % phog_row.num_nonredundant_sequences
  print "Tight ortholog queries: "
  print_query_range(0, cum_queries_for_tight_orthologs)
  print "Tight PHOGRow queries: "
  print_query_range(cum_queries_for_tight_orthologs,
                    cum_queries_for_tight_phogrows)
  print "Custom ortholog queries: "
  print_query_range(cum_queries_for_tight_phogrows, 
                    cum_queries_for_custom_orthologs)
  print "Custom PHOGRow queries: "
  print_query_range(cum_queries_for_custom_orthologs,
                    cum_queries_for_custom_phogrows)
  print "Loose ortholog queries: "
  print_query_range(cum_queries_for_custom_phogrows,
                    cum_queries_for_loose_orthologs)
  print "Loose PHOGRow queries: "
  print_query_range(cum_queries_for_loose_orthologs,
                    cum_queries_for_loose_phogrows)


if __name__ == '__main__':
  main()
