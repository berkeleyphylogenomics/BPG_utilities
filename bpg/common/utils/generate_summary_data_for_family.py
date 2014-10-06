#!/usr/bin/env python

import sys
import re
import os

from pfacts003.phylofacts.models import *
from pfacts003.phylofacts.models import get_dir_of_family_accession

if len(sys.argv) < 2:
    print "This script requires the family accession as an argument"
    sys.exit()

family_id = None

try:
    family_id = int(sys.argv[1][3:])
except TypeError:
    print '''The input must be a family accession like bpg0100039'''
    sys.exit()

family = Family.objects.get(id=family_id)

summary_data = family.get_summary_data()

def get_sequences(family):
    leaves = family.canonical_root_node().get_included_leaves()
    return [leaf.sequence_header.uniprot.uniprot_identifier for leaf in leaves if leaf.sequence_header and
            leaf.sequence_header.uniprot]

def get_consensus_sequence(family_accession):
    consensus_sequence_data = []
    if os.path.exists(get_dir_of_family_accession(family_accession)):
        file_path = os.path.join(get_dir_of_family_accession(family_accession), 
                           '%s.con.fa' % (family_accession))
        if os.path.exists(file_path):
            f = open(file_path)
            f.readline()
            consensus_header ='lcl|%s_consensus' % family_accession
            consensus_sequence = ''.join(f.read().split('\n'))
            consensus_sequence_data = [consensus_header, consensus_sequence]
    return consensus_sequence_data
        
def get_go_data(object):
    go_types = ['BP', 'MF', 'CC']
    for i in range(len(go_types)):
        for go_item in summary_data['go_data'][i][1]:
            output.append("%s\t%s, %s, %s" % (go_types[i], 
                go_item[1], go_item[0], go_item[2]))

def get_ec_data(output):
    ecs = summary_data['ec_links']
    ec_pattern = re.compile(r"""<a href='[^']+'>(.*)</a>""")
    
    for i in range(len(ecs)):
        if ecs[i] != '':
          ec_accession_results = ec_pattern.match(ecs[i][0])
          if ec_accession_results is not None:
              accession = ec_accession_results.group(1)

          # get description
          ec_description_results = ec_pattern.match(ecs[i][1]) 
          if ec_description_results is not None:
              description = ec_description_results.group(1)

          if accession and description:
              output.append("EC\t%s, %s" % (accession.strip(),
                  description.strip()))

def get_pfam_description(output):
    if family.family_type.id == 'C':
        pfam =  family.get_pfams().pop()[0]
        pfam_description = "%s (%s)" % (pfam.name, pfam.accession)
        output.append("PF\t%s" % pfam_description)

def append_sequence_data(header, sequence, entry_code, output):
    start = 0
    chars_per_line = 60
    sequence_length = len(sequence)
    output.append('%s\t>%s' % (entry_code, header))
    while start < sequence_length:
        output.append("\t%s" %  sequence[start:start + chars_per_line])
        start += chars_per_line

# -----------------------------------------------------------------------------
# Start generating output
output = []

# AC (Accession)
output.append("AC\t%s" % family.get_accession())

# CM (Homology Clustering Mode)
clustering_mode = ''
if family.family_type.id == 'C':
    clustering_mode = 'global-local'
else:
    clustering_mode = 'global-global'
output.append("CM\t%s" % clustering_mode) 

# BP (GO Biological Process)
# MF (GO Molecular Function)
# CC (GO Cellular Component)
get_go_data(output)

get_ec_data(output)

get_pfam_description(output)

# NS (Number of Sequences)
output.append("NS\t%d" % len(family.canonical_root_node().get_included_leaves()))

# UD (UniProt Description)
output.append("UD\t%s" % summary_data['description'] )

# AL (Alignment Length)
output.append("AL\t%s" % family.get_alignment_length())

# TD (Taxonomic Distribution)
output.append("FT\t%s" % family.get_taxon())

# SS (Seed Sequence)
try:
  append_sequence_data(family.seed_sequence_header.header,family.seed_sequence_header.sequence.chars, 'SS', output)
except AttributeError:
  output.append("SS\tNot Available")

# CS (Consensus Sequence)
consensus_sequence_header, consensus_sequence = get_consensus_sequence(family.get_accession())
append_sequence_data(consensus_sequence_header, consensus_sequence, 'CS', output)

# MS (Member Sequences):
sequences_per_line = 6
start = 0
sequences = list(set(get_sequences(family)))
sequences.sort()
output.append("MS\t%s" % ', '.join(sequences[start:sequences_per_line]))
start = start + sequences_per_line
while start < len(sequences):
    output.append("\t%s" %  ', '.join(sequences[start:start + sequences_per_line]))
    start = start + sequences_per_line

outf = open(os.path.join(get_dir_of_family_accession(family.get_accession()),
              '%s_summary.dat' % family.get_accession()), 'w')
for line in output:
    outf.write("%s\n" % line)

outf.close()
