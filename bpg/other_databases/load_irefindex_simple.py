#!/usr/bin/env python

from pfacts003.phylofacts.models import UniProt, PPI_DetectionMethod, \
  PPI_SourceDatabase, IRefIndexSimple
import os, re

def get_uniprot_from_alternates(field):
  alts = field.split('|')
  for alt in alts:
    alt_components = alt.split(':')
    db = alt_components[0]
    if db == 'uniprotkb':
      acc = alt_components[1]
      try:
        uniprot = UniProt.objects.get(accession = acc)
      except UniProt.DoesNotExist:
        return None
      return uniprot
  return None

def get_uniprot_partners(fields):
  dbA, accA = fields[0].split(':')
  dbB, accB = fields[1].split(':')
  if dbA == 'uniprotkb':
    try:
      uniprotA = UniProt.objects.get(accession = accA)
    except UniProt.DoesNotExist:
      return None
  else:
    uniprotA = get_uniprot_from_alternates(fields[2])
  if not uniprotA:
    return None
  if dbB == 'uniprotkb':
    try:
      uniprotB = UniProt.objects.get(accession = accB)
    except UniProt.DoesNotExist:
      return None
  else:
    uniprotB = get_uniprot_from_alternates(fields[3])
  if not uniprotB:
    return None
  return (uniprotA, uniprotB)
      
def main():
  # This file should have been made from the file downloaded from iRefIndex with
  # a command like:
  # cat All.mitab.01192011.txt | cut -f 1,2,3,4,7,13,14,53,54 \
  # > All.mitab.cols_1_2_3_4_7_13_14_53_54.txt
  filename = os.path.join('/clusterfs/ohana/external/iRefIndex/current/',
                    'All.mitab.cols_1_2_3_4_7_13_14_53_54.txt')
  mi_re = re.compile('MI:(\d+)\(([^\)]*)\)')
  f = open(filename)
  for line in f.readlines()[1:]:
    fields = line.strip().split('\t')
    partners = get_uniprot_partners(fields)
    if partners:
      uniprotA, uniprotB = partners
      m = mi_re.match(fields[4]) # column 7 of the original file
      method = None
      if m:
        if m.group(1) == '0000':
          methods = PPI_DetectionMethod.objects.filter(short_label = m.group(2))
          if methods:
            method = methods[0]
          else:
            method = PPI_DetectionMethod.objects.create(short_label = m.group(2))
        else:
          methods = PPI_DetectionMethod.objects.filter(
                          mi_ontology_id = m.group(1), short_label = m.group(2))
          if methods:
            method = methods[0]
          else:
            method = PPI_DetectionMethod.objects.create(
                          mi_ontology_id = m.group(1), short_label = m.group(2))
      m = mi_re.match(fields[5]) # column 13 of the original file
      if not m:
        print "Malformed source", line
        continue
      if m.group(1) == '0000':
        sources = PPI_SourceDatabase.objects.filter(source_name = m.group(2))
        if sources:
          source = sources[0]
        else:
          source = PPI_SourceDatabase.objects.create(source_name = m.group(2))
      else:
        sources = PPI_SourceDatabase.objects.filter(mi_ontology_id = m.group(1),
                                                source_name = m.group(2))
        if sources:
          source = sources[0]
        else:
          source = PPI_SourceDatabase.objects.create(mi_ontology_id = m.group(1),
                                                source_name = m.group(2))
      source_id = fields[6].split('|')[0].split(':')[1] # fields[6] is col 14
      interaction_type = fields[7] # column 53 of the original file
      num_participants = fields[8] # column 54 of the original file
      if method:
        IRefIndexSimple.objects.create(uniprot_a = uniprotA, 
                              uniprot_b = uniprotB,
                              ppi_detection_method = method,
                              ppi_source_database = source,
                              interaction_identifier_in_source_db = source_id,
                              interaction_type = interaction_type,
                              num_participants = num_participants)
      else:
        IRefIndexSimple.objects.create(uniprot_a = uniprotA, 
                              uniprot_b = uniprotB,
                              ppi_source_database = source,
                              interaction_identifier_in_source_db = source_id,
                              interaction_type = interaction_type,
                              num_participants = num_participants)
  f.close()

if __name__ == '__main__':
  main()
