#!/usr/bin/env python 

import os.path
import re
import sys
import subprocess

from pfacts003.phylofacts.models import *
from bpg.common.utils.dir_of_family import get_dir_of_family_accession as get_dir

if __name__ == "__main__":

     if len(sys.argv) < 2:
         print "You need to specify a family accession as an arg"
         sys.exit()
 
     family_accession = sys.argv[1]
     try:
         family_id = int(family_accession[3:])
     except ValueError:
         print "You need to specify a family accession as an arg"
         sys.exit()
 
     try:
         family = Family.objects.get(id = family_id)    
         if family.status == "bad":
             raise Family.DoesNotExist
         if family.family_type.id != u'C':
             print "Not a PFAM family"
             sys.exit()
     except Family.DoesNotExist:
         print "No family found for the accession provided"
         sys.exit()
     
     family_dir = get_dir(family_accession)
     realpath = os.path.realpath(family_dir)
     print "realpath = %s" % realpath
     pfam = None
     pfam_accession_pattern = re.compile(r'(PF\d+)\.')
     
     results = pfam_accession_pattern.search(realpath)
     if results is not None:
         pfam_accession = results.group(1)
         print "pfam_accession = %s" % pfam_accession
         try:
             # Pfam version 23 was used for the initial construction of
             # the large, untargeted Pfam library.
             pfam = Pfam.objects.get(accession = pfam_accession,
                 overall_pfam_version = 23)
         except Pfam.DoesNotExist:
             print "Pfam accession %s is missing from the Pfam table" % pfam_accession
         # This is the version of the particular Pfam domain.
         print "pfam.version = %s" % pfam.version
 
     # Case 2: HPylori. 
     hpylori_pattern = re.compile(r'Hpylori26695')
     results = hpylori_pattern.search(realpath)
     if results is not None:
         sub_sequence = os.path.split(realpath)[1]
         print "sub_sequence = %s" % sub_sequence
         gene, start, end = sub_sequence.split("_")
         print [gene, start, end]
         p = subprocess.Popen(['grep', 'gnl|%s' % gene,
             '/clusterfs/ohana/bpg/Hpylori26695/Pfams/domains'],
             stdout = subprocess.PIPE, stderr = subprocess.STDOUT)            
         
         output, error = p.communicate(input=None)
         print output
         lines = output.split('\n')
         domain_name = None
         for line in lines:
             tokens = line.split(',')
             line_start, line_end = tokens[8:10]
             if line_start == start and line_end == end:
                 domain_name = tokens[1]
                 break
         if domain_name:
             print "domain_name: %s" % domain_name
             try:
                 pfam = Pfam.objects.get(name = domain_name,
                     overall_pfam_version = 24)
                 print "pfam.version = %s" % pfam.version 
             except Pfam.DoesNotExist:
                 print "Didn't find Pfam object with the domain name."
                 sys.exit()
         else:
             print "No domain was found in the H. pylori domains file"
             sys.exit() 

     # case 3: YQ's Prokaryote families
     prokaryote_patterns = {
        'ecoli' : (re.compile(r'Ecoli'), '/clusterfs/ohana/bpg/Ecoli_Pfam_books/Ecoli_pfam_subsequences.fa'),
        'mtb'   : (re.compile(r'MTB'),'/clusterfs/ohana/bpg/MTB_Pfam_books/MTB_pfam_subsequences.fa'),
        'qfo'   : (re.compile(r'QuestForOrthologs_eight'), '/clusterfs/ohana/bpg/QuestForOrthologs_eight_prokaryotes/Initial_eight_concatenated_pfam_subsequences.fa')
     }

     for key in prokaryote_patterns:
         results = prokaryote_patterns[key][0].search(realpath)
         if results is not None:
             break
     if results is not None:
         sub_sequence = os.path.split(realpath)[1]
         print "sub_sequence = %s" % sub_sequence
        
         p = subprocess.Popen(['grep', sub_sequence,
            prokaryote_patterns[key][1]],
            stdout = subprocess.PIPE, stderr = subprocess.STDOUT)            
         output, error = p.communicate(input=None)
         print output
         domain_name = output.split(" ,")[1]
         print "pfam_domain = %s" % domain_name
         try:
             pfam = Pfam.objects.get(name = domain_name,
                overall_pfam_version = 24)
             print "pfam.version = %s" % pfam.version 
         except Pfam.DoesNotExist:
             print "Didn't find Pfam object with the domain name."
             sys.exit()

     treefam_pattern = re.compile(r'TreeFamPHOG2')
     result = treefam_pattern.search(realpath)
     if result is not None:
        domain_name = os.path.split(realpath)[1]
        print "pfam_domain = %s" % domain_name
        try:
          pfam = Pfam.objects.get(name = domain_name, 
                                  overall_pfam_version = 24)
          print "pfam.version = %s" % pfam.version
        except Pfam.DoesNotExist:
          print "Didn't find Pfam object with the domain name."
          sys.exit()

     if pfam is None:
        print "%s FAIL: no Pfam object associated to conserved region family" \
              % family_accession
     else:
        canonical_root_consensus = TreeNodeConsensus.objects.get(
                                    tree_node = family.canonical_root_node(), 
                                    method = 'hmm')
        canonical_root_consensus.hmm_consensus.hmm.pfam = pfam
        canonical_root_consensus.hmm_consensus.hmm.save()
        print "%s SUCCESS: Saved association to pfam %s %s to the database" \
              % (family_accession, pfam.accession, pfam.name)
