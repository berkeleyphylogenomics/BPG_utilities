#!/usr/bin/env python

import sys
from pfacts003.phylofacts.models import *
from optparse import OptionParser
from pfacts003.utils.id_patterns import is_uniprot_accession_format
from pfacts003.utils.json_data_objects import _annotate_nhx_element
from Bio import SwissProt

def indented_write(outf, indentation_level, buffer):
  for i in range(0, indentation_level):
    outf.write('    ')
  outf.write(buffer)

def write_swissprot_annotations(outf, indentation_level, uniprot, uniprot_f):
  uniprot_dat_indices = UniProtDatIndex.objects.filter(uniprot = uniprot)
  for uniprot_dat_index in uniprot_dat_indices:
    if uniprot_dat_index.uniprot_accession == uniprot.accession:
      break
  uniprot_f.seek(uniprot_dat_index.file_char)
  record = SwissProt.parse(uniprot_f).next()
  indented_write(outf, indentation_level+1, 
                  "Length: %d\n" % record.sequence_length)
  if len(record.gene_name) > 0:
    for name_spec in record.gene_name.replace('\n',' ').split('; '):
      name_type, names = name_spec.split('=')
      indented_write(outf, indentation_level+1, 
                      "%s: %s\n" % (name_type, names))
  for keyword in record.keywords:
    indented_write(outf, indentation_level+1, 'Keyword: %s\n' % keyword)
  for comment in record.comments:
    if comment[0:5] == '-----':
      continue
    components = comment.replace(':\n',': ').split(': ')
    comment_type = components[0]
    comment_lines = ': '.join(components[1:]).split('\n')
    indented_write(outf, indentation_level+1, "%s:\n" % comment_type)
    for line in comment_lines:
      indented_write(outf, indentation_level+2, "%s\n" % line)
  for cross_reference in record.cross_references:
    indented_write(outf, indentation_level+1, 
                  "%s: %s\n" % (cross_reference[0],
                                '; '.join(cross_reference[1:])))

def write_annotations(outf, indentation_level, sequence_header, annotations,
                      uniprot_f, kegg_maps_of_ec):
  indented_write(outf, indentation_level, 
                  "%s {\n" % sequence_header.identifier())
  inSwissProt = False
  for annotation_type in annotations[str(sequence_header.id)]:
    if type(annotations[str(sequence_header.id)][annotation_type]) == str or \
        type(annotations[str(sequence_header.id)][annotation_type]) == unicode:
      output = annotations[str(sequence_header.id)][annotation_type]
    elif annotation_type == 'in_swissprot_f':
     inSwissProt = annotations[str(sequence_header.id)][annotation_type] != 0
     if inSwissProt:
       output = str(inSwissProt) 
     else:
       output = 'N/A'
    elif annotation_type == 'ncbi_taxid':
      output = "%d" % annotations[str(sequence_header.id)][annotation_type]
    else: 
      output = ';'.join([str(annotation) for annotation in annotations[
                                    str(sequence_header.id)][annotation_type]])
    if output != 'N/A':
      indented_write(outf, indentation_level+1, "%s: %s\n" % (annotation_type,
                                                              output))
  if inSwissProt:
    write_swissprot_annotations(outf, indentation_level,
                                sequence_header.uniprot, uniprot_f)
  if sequence_header.uniprot:
    pfam_hits = sequence_header.uniprot.get_pfam_hits()
    for pfam_hit in pfam_hits:
      indented_write(outf, indentation_level+1,
      "Pfam: Domain %s Residues %s-%s E-value %s\n"
      % (pfam_hit.pfamA.description, pfam_hit.ali_start, pfam_hit.ali_end,
      pfam_hit.sequence_evalue_score))
    partner_uniprots = sequence_header.uniprot.get_interacting_partners()
    if len(partner_uniprots) > 0:
      indented_write(outf, indentation_level+1, 
          "Interacting partners: %s\n" % ','.join([uniprot.accession for uniprot
                                                in partner_uniprots]))
  kegg_maps = set()
  if 'ec' in annotations[str(sequence_header.id)]:
    ec_strs = annotations[str(sequence_header.id)]['ec'].split(';')
    for ec_str in ec_strs:
      if ec_str in kegg_maps_of_ec:
        kegg_maps = kegg_maps | kegg_maps_of_ec[ec_str]

  for kegg_map in kegg_maps:
    indented_write(outf, indentation_level+1, "KEGG %s\n" % kegg_map.__str__())

  indented_write(outf, indentation_level,
                  "} %s\n" % sequence_header.identifier())

def main():
  usage = "%prog [options] family_accession"
  opt_parser = OptionParser(usage=usage)
  (options, args) = opt_parser.parse_args()
  if len(args) != 1:
    opt_parser.error('Incorrect number of arguments')
  if len(args[0]) != 10 or args[0][0:3] != 'bpg':
    opt_parser.error('Family accession is not a valid bpg accession.')
  try:
    family_id = int(args[0][3:])
  except ValueError:
    opt_parser.error('Family accession is not a valid bpg accession.')
  families = Family.objects.filter(id = family_id).exclude(status='bad')
  if families:
    family = families[0]
  else:
    opt_parser.error('Family accession is not a valid bpg accession.')
  trees = Tree.objects.filter(family = family, method = 'ml')
    
  kegg_map_ecs = KEGG_Map_EC.objects.all()
  kegg_maps_of_ec = {}
  for kegg_map_ec in kegg_map_ecs:
    ec_str = kegg_map_ec.ec.__unicode__()
    if ec_str not in kegg_maps_of_ec:
      kegg_maps_of_ec[ec_str] = set()
    kegg_maps_of_ec[ec_str].add(kegg_map_ec.kegg_map)

  uniprot_f = open('/clusterfs/ohana/external/UniProt/current/uniprot.dat')
  for tree in trees:
    max_right_id = 0
    leaves = TreeNode.objects.filter(tree = tree, 
                          sequence_header__isnull = False).order_by('left_id')
    for leaf in leaves:
      max_right_id = max(max_right_id, leaf.left_id)
      
    annotations = _annotate_nhx_element(tree.family.get_accession(),
                                  [leaf.sequence_header.id for leaf in leaves])
    left_phogs = TreeNode.objects.filter(tree = tree, 
                    duplication_distance__isnull = False).order_by('left_id')
    right_phogs = TreeNode.objects.filter(tree = tree, 
                    duplication_distance__isnull = False).order_by('right_id')
    sys.stdout.write("%s {\n" % tree.family.get_accession())
    indentation_level = 1
    left_id = 0
    leaf_idx = 0
    left_phog_idx = 0
    right_phog_idx = 0
    while leaf_idx < len(leaves) or left_phog_idx < len(left_phogs) \
        or right_phog_idx < len(right_phogs):
      left_id = -1
      next_leaf = None
      next_left_phog = None
      next_right_phog = None
      if leaf_idx < len(leaves):
        next_leaf = leaves[leaf_idx]
        left_id = next_leaf.left_id
      if left_phog_idx < len(left_phogs):
        next_left_phog = left_phogs[left_phog_idx]
        if left_id >= 0:
          left_id = min(left_id, next_left_phog.left_id)
        else:
          left_id = next_left_phog.left_id
      if right_phog_idx < len(right_phogs):
        next_right_phog = right_phogs[right_phog_idx]
        if left_id >= 0:
          left_id = min(left_id, next_right_phog.right_id)
        else:
          left_id = next_right_phog.right_id
      if next_leaf and left_id == next_leaf.left_id:
        leaf_idx += 1
        write_annotations(sys.stdout, indentation_level,
                          next_leaf.sequence_header, annotations, uniprot_f,
                          kegg_maps_of_ec)
      elif next_left_phog and left_id == next_left_phog.left_id:
        if next_left_phog.duplication_distance > \
            next_left_phog.greatest_duplication_distance_of_maximal_descendant:
          if next_left_phog.greatest_duplication_distance_of_maximal_descendant <= 0.0:
            indented_write(sys.stdout, indentation_level,
              "%s Thresholds [0.0,%g] {\n" % (next_left_phog.get_accession(),
              next_left_phog.duplication_distance))
          else:
            indented_write(sys.stdout, indentation_level,
              "%s Thresholds (%g,%g] {\n" % (next_left_phog.get_accession(),
              next_left_phog.greatest_duplication_distance_of_maximal_descendant,
              next_left_phog.duplication_distance))
          indentation_level += 1
        left_phog_idx += 1
      elif next_right_phog and left_id == next_right_phog.right_id:
        if next_right_phog.duplication_distance > \
            next_right_phog.greatest_duplication_distance_of_maximal_descendant:
          indentation_level -= 1
          indented_write(sys.stdout, indentation_level,
            "} %s\n" % next_right_phog.get_accession())
        right_phog_idx += 1
    sys.stdout.write("} %s\n" % tree.family.get_accession())
  uniprot_f.close()

if __name__ == '__main__':
  main()
