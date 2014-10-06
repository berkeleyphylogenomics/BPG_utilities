#!/usr/bin/env python

import os, re
from optparse import OptionParser
from pfacts003.phylofacts.models import Pfam, HMM

def main():
  parser = OptionParser(usage='%prog')
  (options, args) = parser.parse_args()

  overall_pfam_version = 23.0
  hmm_type = 'long'
  f = open("/clusterfs/ohana/external/pfam/Pfam23.0/Pfam_ls")
  for line in f.readlines():
    line_type = line[0:5].strip()
    if line_type == 'NAME':
      name = line[6:].strip()
    elif line_type == 'ACC':
      accession, version_str = line[6:].strip().split('.')
      version = int(version_str)
    elif line_type == 'DESC':
      description = line[6:].strip()
    elif line_type == 'LENG':
      hmmlength = int(line[6:].strip())
    elif line_type == '//':
      pfam_objects = Pfam.objects.filter(accession__exact = accession,
                            overall_pfam_version__exact = overall_pfam_version)
      if pfam_objects:
        pfam = pfam_objects[0]
      else:
        pfam = Pfam.objects.create(accession = accession,
                                  name = name,
                                  description = description,
                                  hmmlength = hmmlength,
                                  hmm_type = hmm_type,
                                  version = version,
                                  overall_pfam_version = overall_pfam_version)
      hmm_objects = HMM.objects.filter(pfam__exact = pfam, 
                                        hmm_type__exact = 'HMMER2',
                                        method__exact = '3rd_party')
      if not hmm_objects:
        hmm = HMM.objects.create(length = hmmlength,
                                  hmm_type = 'HMMER2',
                                  method = '3rd_party',
                                  pfam = pfam)
      accession = None
      name = None
      description = None
      hmmlength = None
      version = None

if __name__ == '__main__':
  main()
