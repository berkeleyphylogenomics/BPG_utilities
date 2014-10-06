#!/usr/bin/env python

import django
from pfacts003.phylofacts import models
import sys

global verbose

def update_left_right_ids(taxon, left_id, level):
  """Execute a modified pre-order tree traversal and store the left & right ids.

  See for instance 
  http://articles.sitepoint.com/article/hierarchical-data-database/2
  or
  http://www.oreillynet.com/pub/a/network/2002/1/27/bioconf.html?page=1
  """
  global verbose
  if verbose:
    for i in xrange(level):
      print " ",
    print "Taxon: id: %06d  left_id: %d" % (taxon.id, left_id)
  sys.stdout.flush()
  taxon.left_id = left_id
  taxon.save()
  right_id = left_id + 1
  for child in taxon.children.all():
    right_id = update_left_right_ids(child, right_id, level + 1)
  if verbose:
    for i in xrange(level):
      print " ",
    print "Taxon: id: %06d right_id: %d" % (taxon.id, right_id)
  sys.stdout.flush()
  taxon.right_id = right_id
  taxon.save()
  return (right_id + 1)
  
  
def main():
  global verbose
  verbose = False
  root_taxon = models.UniProtTaxonomy.objects.get(id__exact = 1)
  right_id = update_left_right_ids(root_taxon, 1, 0)

if __name__=="__main__":
  main()
