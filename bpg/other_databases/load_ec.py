#!/usr/bin/env python

from optparse import OptionParser
from pfacts003.phylofacts.models import EC

def main():
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()

  f = open("/clusterfs/ohana/external/IntEnz/current/enzclass.txt")
  lines = f.readlines()
  f.close()

  for line in lines:
    description = line[11:].strip()
    enzyme_classification = line[0:11]
    fields = [field.strip() for field in enzyme_classification.split('.')]
    # Skip the initial lines describing the file
    if len(fields) < 4:
      continue
    class_number = int(fields[0])
    if fields[1] == '-':
      ec_objects = EC.objects.filter(class_number__exact = class_number,
                                      subclass_number__isnull = True,
                                      subsubclass_number__isnull = True,
                                      enzyme_number__isnull = True)
      if ec_objects:
        ec = ec_objects[0]
        ec.description = description
        ec.save()
      else:
        ec = EC.objects.create(class_number = class_number, 
                                description = description)
    else:
      subclass_number = int(fields[1])
      if fields[2] == '-':
        ec_objects = EC.objects.filter(class_number__exact = class_number,
                                      subclass_number__exact = subclass_number,
                                      subsubclass_number__isnull = True,
                                      enzyme_number__isnull = True)
        if ec_objects:
          ec = ec_objects[0]
          ec.description = description
          ec.save()
        else:
          ec = EC.objects.create(class_number = class_number, 
                                  subclass_number = subclass_number,
                                  description = description)
      else:
        subsubclass_number = int(fields[2])
        if fields[3] == '-':
          ec_objects = EC.objects.filter(class_number__exact = class_number,
                                subclass_number__exact = subclass_number,
                                subsubclass_number__exact = subsubclass_number,
                                enzyme_number__isnull = True)
          if ec_objects:
            ec = ec_objects[0]
            ec.description = description
            ec.save()
          else:
            ec = EC.objects.create(class_number = class_number, 
                                    subclass_number = subclass_number,
                                    subsubclass_number = subsubclass_number,
                                    description = description)
        else:
          # Actually, we don't expect this ever to happen in enzclass.txt
          enzyme_number = int(fields[3])
          ec_objects = EC.objects.filter(class_number__exact = class_number,
                                subclass_number__exact = subclass_number,
                                subsubclass_number__exact = subsubclass_number,
                                enzyme_number__exact = enzyme_number)
          if ec_objects:
            ec = ec_objects[0]
            ec.description = description
            ec.save()
          else:
            ec = EC.objects.create(class_number = class_number, 
                                    subclass_number = subclass_number,
                                    subsubclass_number = subsubclass_number,
                                    enzyme_number = enzyme_number,
                                    description = description)

  f = open("/clusterfs/ohana/external/IntEnz/current/enzyme.dat")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line[0:2] == 'ID':
      class_number, subclass_number, subsubclass_number, enzyme_number \
        = [int(x) for x in line[2:].strip().split('.')]
    elif line[0:2] == 'DE':
      description = line[2:].strip()
      ec_objects = EC.objects.filter(class_number__exact = class_number,
                            subclass_number__exact = subclass_number,
                            subsubclass_number__exact = subsubclass_number,
                            enzyme_number__exact = enzyme_number)
      if ec_objects:
        ec = ec_objects[0]
        ec.description = description
        ec.save()
      else:
        ec = EC.objects.create(class_number = class_number, 
                                subclass_number = subclass_number,
                                subsubclass_number = subsubclass_number,
                                enzyme_number = enzyme_number,
                                description = description)
  

if __name__ == '__main__':
  main()
