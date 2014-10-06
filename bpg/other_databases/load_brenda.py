#!/usr/bin/env python

import re
from optparse import OptionParser
from pfacts003.phylofacts.models import UniProtDatIndex,UniProt,UniProtEC,EC

def main():
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()

  f = open("/clusterfs/ohana/external/BRENDA/current/brenda.txt")
  lines = f.readlines()
  f.close()

  ec_number = ""
  ec_object = None
  for line in lines:
    if line[0:3] == 'ID\t':
      ec_number = line[2:].strip().split()[0]
      fields = [field.strip() for field in ec_number.split('.')]
      class_number = int(fields[0])
      subclass_number = int(fields[1])
      subsubclass_number = int(fields[2])
      try:
        enzyme_number = int(fields[3])
      except ValueError:
        # This is probably a preliminary BRENDA-supplied EC number.
        # We have "preliminary" numbers in the EC table, but these come from
        # UniProt.  There is no guarantee that the preliminary BRENDA-supplied
        # EC numbers match with the UniProt ones, so we will just skip these
        # cases.
        ec_object = None
        continue
      ec_objects = EC.objects.filter(class_number__exact = class_number,
                                subclass_number__exact = subclass_number,
                                subsubclass_number__exact = subsubclass_number,
                                enzyme_number__exact = enzyme_number)
      if ec_objects:
        ec_object = ec_objects[0]
      else:
        ec_object = None
        print "EC %s missing from database" % ec_number
    elif line[0:3] == 'PR\t' and ec_object:
      words = line[2:].strip().split()
      if words[-2] == 'UniProt' or words[-2] == 'SwissProt' \
          or words[-2] == 'TrEMBL':
        uniprot = None
        uniprot_accession = words[-3]
        try:
          uniprot = UniProt.objects.get(accession__exact = uniprot_accession)
        except UniProt.DoesNotExist:
          try:
            uniprot_dat_index = UniProtDatIndex.objects.get(
                                  uniprot_accession__exact = uniprot_accession)
            uniprot = uniprot_dat_index.uniprot
          except UniProtDatIndex.DoesNotExist:
            pass
        if uniprot:
          uniprot_ec_objects = UniProtEC.objects.filter(
                                  uniprot__exact = uniprot,
                                  ec__exact = ec_object)
          if uniprot_ec_objects:
            uniprot_ec_object = uniprot_ec_objects[0]
          else:
            uniprot_ec_object = UniProtEC.objects.create(uniprot = uniprot,
                                                        ec = ec_object)
          uniprot_ec_object.is_in_brenda_f = True
          uniprot_ec_object.save()

if __name__ == '__main__':
  main()
