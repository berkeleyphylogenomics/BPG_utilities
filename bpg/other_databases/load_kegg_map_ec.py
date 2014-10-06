#!/usr/bin/env python

from optparse import OptionParser
from pfacts003.phylofacts.models import KEGG_Map, KEGG_Map_EC, EC

kegg_map_object_of_id = {}

def get_kegg_map_object(map_id):
  if map_id not in kegg_map_object_of_id:
    kegg_map_object, created = KEGG_Map.objects.get_or_create(id=int(map_id))
    kegg_map_object_of_id[map_id] = kegg_map_object
  return kegg_map_object_of_id[map_id]

def main():
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()

  f = open("/clusterfs/ohana/external/KEGG/current/ec_map.tab")
  lines = f.readlines()
  f.close()
  for line in lines:
    ec_number, map_line = line.strip().split('\t')
    ec_objects = None
    maps = map_line.split()
    fields = ec_number.split('.')
    try:
      class_number = int(fields[0])
      if len(fields) == 1 or fields[1] == '-':
        ec_objects = EC.objects.filter(class_number__exact = class_number,
                                      subclass_number__isnull = True,
                                      subsubclass_number__isnull = True,
                                      enzyme_number__isnull = True)
      else:
        subclass_number = int(fields[1])
        if len(fields) == 2 or fields[2] == '-':
          ec_objects = EC.objects.filter(class_number__exact = class_number,
                                      subclass_number__exact = subclass_number,
                                      subsubclass_number__isnull = True,
                                      enzyme_number__isnull = True)
        else:
          subsubclass_number = int(fields[2])
          if len(fields) == 3 or fields[3] == '-':
            ec_objects = EC.objects.filter(class_number__exact = class_number,
                                subclass_number__exact = subclass_number,
                                subsubclass_number__exact = subsubclass_number,
                                enzyme_number__isnull = True)
          else:
            enzyme_number = int(fields[3])
            ec_objects = EC.objects.filter(class_number__exact = class_number,
                                subclass_number__exact = subclass_number,
                                subsubclass_number__exact = subsubclass_number,
                                enzyme_number__exact = enzyme_number,
                                is_preliminary_f__exact = False)
    except ValueError:
      # There are various bugs in the KEGG ec_map.tab
      # so there's no way of telling what this is actually supposed to be.
      # Just skip it.
      ec_objects = None

    if ec_objects:
      ec = ec_objects[0]
      for map in maps:
        kegg_map = get_kegg_map_object(map)
        kegg_map_ec_objects = KEGG_Map_EC.objects.filter(ec__exact = ec,
                                                    kegg_map__exact = kegg_map)
        if not kegg_map_ec_objects:
          KEGG_Map_EC.objects.create(ec = ec, kegg_map = kegg_map)

if __name__ == '__main__':
  main()
