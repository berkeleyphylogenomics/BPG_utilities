#!/usr/bin/env python

from pfacts003.phylofacts.models import SequenceHeader, UniProt, Family

def main():
  f = open('/home/yshen/seed_for_still_failed_good_families_20110808.txt')
  lines = [line.strip() for line in f.readlines()]
  f.close()

  family_accession = None
  header = ''
  for line in lines:
    if len(line) == 0:
      family_accession = None
      header = ''
    elif family_accession is not None:
      if line != 'No seed' and not done_family:
        if header == '':
          header = line[1:]
        else:
          header = ' '.join([header, line[1:]])
        identifier = header.split()[0]
        accession = identifier.split('|')[1]
        uniprot =  UniProt.objects.get(accession = accession)
        sequence_headers = SequenceHeader.objects.filter(uniprot = uniprot)
        if len(sequence_headers) > 1:
          sequence_headers = sequence_headers.filter(header = header)
        if len(sequence_headers) == 1:
          print "Assigning seed_sequence_header %d for family %s" \
              % (sequence_headers[0].id, family_accession)
          family = Family.objects.get(id = int(family_accession[3:]))
          family.seed_sequence_header = sequence_headers[0]
          family.save()
    else:
      family_accession = line.split(',')[0]
      done_family = False
      header = ''

if __name__ == '__main__':
  main()
