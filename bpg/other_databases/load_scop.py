#!/usr/bin/env python

import os, re
from optparse import OptionParser
from pfacts003.phylofacts.models import SCOP, PDB_SCOP, PDB, PDB_Chain

def main():
  parser = OptionParser(usage='%prog <scop_description_path>')
  (options, args) = parser.parse_args()
  if len(args) < 1:
    parser.error("Path to SCOP description is required.")
  scop_description_path = args[0]
  if not os.path.exists(scop_description_path):
    parser.error("SCOP description file %s not found." % scop_description_path)
  scop_version = float(os.path.split(scop_description_path)[1].split('_')[1])
  f = open(scop_description_path)
  lines = [line.rstrip() for line in f.readlines()]
  f.close()
  range_re = re.compile("(-?[0-9]*)[ABCPS]?-([0-9]*)[ABCPS]?")
  for line in lines:
    if line[0] == '#':
      continue
    fields = line.split('\t')
    scop_specification = fields[2]
    if fields[1] == 'cl':
      class_letter = scop_specification
      scop_objects = SCOP.objects.filter(class_letter__exact = class_letter,
                                        fold_number__isnull = True,
                                        superfamily_number__isnull = True,
                                        family_number__isnull = True)
      if scop_objects:
        scop = scop_objects[0]
      else:
        scop = SCOP.objects.create(class_letter = class_letter)
      description = fields[4]
      scop.description = description
      scop.save()
    elif fields[1] == 'cf':
      class_letter, fold_string = scop_specification.split('.')
      fold_number = int(fold_string)
      scop_objects = SCOP.objects.filter(class_letter__exact = class_letter,
                                        fold_number__exact = fold_number,
                                        superfamily_number__isnull = True,
                                        family_number__isnull = True)
      if scop_objects:
        scop = scop_objects[0]
      else:
        scop = SCOP.objects.create(class_letter = class_letter,
                                    fold_number = fold_number)
      description = fields[4]
      scop.description = description
      scop.save()
    elif fields[1] == 'sf':
      class_letter, fold_string, superfamily_string \
          = scop_specification.split('.')
      fold_number = int(fold_string)
      superfamily_number = int(superfamily_string)
      scop_objects = SCOP.objects.filter(class_letter__exact = class_letter,
                                fold_number__exact = fold_number,
                                superfamily_number__exact = superfamily_number,
                                family_number__isnull = True)
      if scop_objects:
        scop = scop_objects[0]
      else:
        scop = SCOP.objects.create(class_letter = class_letter,
                                    fold_number = fold_number,
                                    superfamily_number = superfamily_number)
      description = fields[4]
      scop.description = description
      scop.save()
    elif fields[1] == 'fa':
      class_letter, fold_string, superfamily_string, family_string \
          = scop_specification.split('.')
      fold_number = int(fold_string)
      superfamily_number = int(superfamily_string)
      family_number = int(family_string)
      scop_objects = SCOP.objects.filter(class_letter__exact = class_letter,
                                fold_number__exact = fold_number,
                                superfamily_number__exact = superfamily_number,
                                family_number__exact = family_number)
      if scop_objects:
        scop = scop_objects[0]
      else:
        scop = SCOP.objects.create(class_letter = class_letter,
                                    fold_number = fold_number,
                                    superfamily_number = superfamily_number,
                                    family_number = family_number)
      description = fields[4]
      scop.description = description
      scop.save()
    elif fields[1] == 'px':
      class_letter, fold_string, superfamily_string, family_string \
          = scop_specification.split('.')
      fold_number = int(fold_string)
      superfamily_number = int(superfamily_string)
      family_number = int(family_string)
      scop_objects = SCOP.objects.filter(class_letter__exact = class_letter,
                                fold_number__exact = fold_number,
                                superfamily_number__exact = superfamily_number,
                                family_number__exact = family_number)
      if scop_objects:
        scop = scop_objects[0]
      else:
        scop = SCOP.objects.create(class_letter = class_letter,
                                    fold_number = fold_number,
                                    superfamily_number = superfamily_number,
                                    family_number = family_number)
      scop_identifier = fields[3]
      pdb_identifier, rest_specification = fields[4].split()
      pdb_objects = PDB.objects.filter(id__exact = pdb_identifier)
      if pdb_objects:
        pdb = pdb_objects[0]
      else:
        # The PDB entry may have been deleted, so skip
        continue
      chain_specifications = rest_specification.split(',')
      # FIXME RSD 2010-02-22: Although a single SCOP identifier can be
      # associated with multiple chains, only the last one is added due to the
      # unique constraint on scop_identifier.  The unique constraint may need
      # to be dropped, but then this has to be done carefully, in such a way
      # that reading the same SCOP file twice does not create duplicate
      # entries.
      for chain_specification in chain_specifications:
        part_specification = chain_specification.split(':')
        chain_id = part_specification[0]
        pdb_chain_objects = PDB_Chain.objects.filter(pdb__exact = pdb,
                                            chain_id__exact = chain_id)
        if pdb_chain_objects:
          pdb_chain = pdb_chain_objects[0]
        else:
          # The chain may have been deleted from the PDB entry, so skip
          continue
        pdb_scop_objects = PDB_SCOP.objects.filter(scop_identifier__exact 
                                                              = scop_identifier)
        if pdb_scop_objects:
          # Uniqueness of scop_identifier guarantees there is only one
          pdb_scop = pdb_scop_objects[0]
          pdb_scop.pdb_chain = pdb_chain
          pdb_scop.scop = scop
          pdb_scop.scop_version = scop_version
        else:
          pdb_scop = PDB_SCOP.objects.create(pdb_chain = pdb_chain,
                                        scop = scop,
                                        scop_identifier = scop_identifier,
                                        scop_version = scop_version)
        if len(part_specification) > 1 and len(part_specification[1]) > 0:
          m = range_re.match(part_specification[1])
          start_residue = int(m.group(1))
          end_residue = int(m.group(2))
          pdb_scop.start_residue = start_residue
          pdb_scop.end_residue = end_residue
        pdb_scop.save()

if __name__ == '__main__':
  main()
