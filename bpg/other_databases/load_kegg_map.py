#!/usr/bin/env python

from optparse import OptionParser
from pfacts003.phylofacts.models import KEGG_Map

def main():
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()

  f = open("/clusterfs/ohana/external/KEGG/current/map_title.tab")
  lines = f.readlines()
  f.close()
  for line in lines:
    map_id, map_title = line.strip().split('\t')
    print map_id, map_title
    kegg_map, created = KEGG_Map.objects.get_or_create(id = int(map_id))
    kegg_map.title = map_title
    kegg_map.save()

if __name__ == '__main__':
  main()
