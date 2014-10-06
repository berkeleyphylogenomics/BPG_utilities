#!/bin/env python
##!/usr/bin/mod9v3   This line doesn't work on ohana.
"""
superpose_model_on_true_structure.py
Written by Ruchira Datta March 2008; modified by T. Farrah July 2008

Given:
The structure of a target has been predicted through
homology modeling using a template. The true structure of the
target is known.

Using Modeller from Sali lab:
Superpose the target model on the true structure of the target.
Compute various measures of the model quality.
"""

ohana = 0
salilab = 1

import os, sys, string, commands
if salilab:
        sys.path = ['/netapp/sali/ruchira/pylib'] + sys.path
from modeller import *
from array import array
from matchmaker.shmm_shmm_lib import *

global _PDB_website
_PDB_website = "www.rcsb.org/pdb/files"

def superpose_and_compute_gdt_ts(mdl, mdl2,  aln):
  """ Compute GDT_TS by colling superpose() in Modeller
  """
  # We only want to look at alpha carbon atoms
  atmsel = selection(mdl).only_atom_types('CA')
  # Superpose four times using different rms cutoffs.
  # Each call should yield the same superposition, but
  # each result will contain a different count of how many atoms
  # were within the cutoff.
  r1 = atmsel.superpose(mdl2, aln, rms_cutoff=1, refine_local=False,
  superpose_refine=False)
  r2 = atmsel.superpose(mdl2, aln, rms_cutoff=2, refine_local=False,
  superpose_refine=False)
  r4 = atmsel.superpose(mdl2, aln, rms_cutoff=4, refine_local=False,
  superpose_refine=False)
  r8 = atmsel.superpose(mdl2, aln, rms_cutoff=8, refine_local=False,
  superpose_refine=False)
  # Compute GDT_TS
  gdt_ts = (r1.num_equiv_cutoff_pos + r2.num_equiv_cutoff_pos +
	    r4.num_equiv_cutoff_pos + r8.num_equiv_cutoff_pos)
  num_residues = len(aln[1].residues)
  gdt_ts = gdt_ts / (4.0 * num_residues)
  # Return one of the superpositions, plus all the various numbers.
  return ( r1,  num_residues, r1.num_equiv_cutoff_pos,
                                r2.num_equiv_cutoff_pos,
                                r4.num_equiv_cutoff_pos, 
                                r8.num_equiv_cutoff_pos, gdt_ts)


def write_ali_file(seed_id, true_struct_seq, model_seq, ali_file,
                    starting_residue, ending_residue):
  """Write a Modeller format alignment file for seed vs. self

  True structure  must be first sequence written to file
  """
  chain_id = seed_id[4:5]
  pdb_id = seed_id[0:4]
  ali_fp = open(ali_file, "w")
  ali_fp.write(">P1;%s\n" % seed_id)
  ali_fp.write("structure:%s:FIRST:%s:LAST:%s:undefined:undefined:-1.00:-1.00\n" % (pdb_id, chain_id, chain_id))
  ali_fp.write(true_struct_seq.tostring() + "*\n")
  ali_fp.write(">P1;%s_model\n" % seed_id)
  ali_fp.write("structureM:%s:%s:%s:%s:%s::: 0.00: 0.00\n" % (pdb_id,
                starting_residue, chain_id, ending_residue, chain_id))
  ali_fp.write(model_seq.tostring() + "*\n")
  ali_fp.close()


def create_target_alignment_with_self_for_modeller(alignment_path,
     seed_id, model_dir, target_id, env, unaligned_true_seq):
  """Create .ali file for target seq aligned with self, read into Modeller

  Each residue is aligned with itself, except for residues
   that were not modeled (i.e. were aligned with gap characters)"""

  true_struct_seq = array('c')
  model_seq = array('c')

  # read Matchmaker alignment file used to model target against template
  """
  # 7/15/08: replace with function that can read compressed alignment files
  afa_fp = open(alignment_path)
  afa_lines = afa_fp.readlines()
  afa_fp.close()
  # figure out which alignment line is target (assume other is template)
  # and copy alignment lines into strings
  if afa_lines[0].rstrip()[1:] == target_id:
    target_line = afa_lines[1].strip().upper()
    template_line = afa_lines[3].strip().upper()
  elif afa_lines[2].rstrip()[1:] == target_id:
    target_line = afa_lines[3].strip().upper()
    template_line = afa_lines[1].strip().upper()
  """
  (seed1, seq1), (seed2, seq2) = read_alignment_file(alignment_path)
  # figure out which alignment line is target (assume other is template)
  # and copy alignment lines into strings
  if seed1 == target_id:
    target_line = seq1.upper()
    template_line = seq2.upper()
  elif seed2 == target_id:
    target_line = seq2.upper()
    template_line = seq1.upper()
  else:
    print "Couldn't find %s in %s." % (target_id, alignment_path)
    sys.exit(0)
  # replace dots with dashes in alignment lines
  target_line = string.replace(target_line, '.', '-')
  template_line = string.replace(template_line, '.', '-')
  unaligned_target_line = string.replace(target_line, '-', '')

  i = 0
  j = 0
  starting_residue = 'FIRST'
  ending_residue = 'LAST'
  while i < len(unaligned_true_seq) and j < len(unaligned_target_line):
    if unaligned_true_seq[i] == unaligned_target_line[j]:
          true_struct_seq.fromstring(unaligned_true_seq[i])
          model_seq.fromstring(unaligned_target_line[j])
          i = i + 1
          j = j + 1
    else:
          if len(unaligned_true_seq) > len(unaligned_target_line):
                true_struct_seq.fromstring(unaligned_true_seq[i])
                model_seq.fromstring('-')
                i = i + 1
          else:
                true_struct_seq.fromstring('-')
                model_seq.fromstring(unaligned_target_line[j])
                j = j + 1
                starting_residue = str(2)

  if i < len(unaligned_true_seq):
    ending_residue = str(i)
  while i < len(unaligned_true_seq):
          true_struct_seq.fromstring(unaligned_true_seq[i])
          model_seq.fromstring('-')
          i = i + 1
  while j < len(unaligned_target_line):
          true_struct_seq.fromstring('-')
          model_seq.fromstring(unaligned_target_line[j])
          j = j + 1

  # write .ali file using self-self alignment constructed just above.
  alignment_basename = os.path.splitext(os.path.split(alignment_path)[1])[0]
  ali_file = os.path.join(model_dir, "%s_%s.ali" % (target_id, seed_id))
  write_ali_file(seed_id, true_struct_seq, model_seq, ali_file,
                starting_residue, ending_residue)
  # pass this alignment filename into the modelling environment
  aln = alignment(env, file=(ali_file), 
		  align_codes = ( seed_id, "%s_model" % seed_id))
  return aln


def main():
  # process command line and initialize stuff

  # Assume we are running this in the same directory as
  # the model (.pdb), the Matchmaker alignment (.afa),
  # and the Modeller alignment (.ali)
  # Optionally, <seed_id>.pdb (true structure) exists in same dir,
  #  else it will be obtained from RCSB.
  # We need to provide the PDB id for the true structure and the
  #  name of the sequence when it was modelled.
  # PDB file for true structure must contain exactly the same residues
  #  as the model.

  if len(sys.argv) < 3:
    print "Usage: %s <pdb_id_of_true_structure> <target_id> " % sys.argv[0] \
                   + " [<model_file>] [<alignment_file>]"
    sys.exit(0)

  seed_id = sys.argv[1]
  target_id = sys.argv[2]

  if len(sys.argv) > 3:
    model_file = sys.argv[3]
  else: model_file = glob.glob("*.B99990001.pdb")[0]
  model_dir = os.path.split(model_file)[0]
  # alignment_path file should fasta format file with exactly 4 lines
  # as produced by Matchmaker
  if len(sys.argv) > 4:
    alignment_path = sys.argv[4]
  else: alignment_path = glob.glob("*.afa")[0]

  # create a new modeler environment
  env = environ()

  # get PDB file from RCSB if it isn't already there
  global _PDB_website
  pdb_filename = os.path.join(single_seqs_dir(), target_id, "%s.pdb" %
                              seed_id[0:4])
  if not os.path.exists(pdb_filename):
    if ohana:
      original_pdb_filename \
          = '/clusterfs/ohana/external/pdb_SHMM-SHMM/pdb%s.ent' % seed_id[0:4]
    else:
      original_pdb_filename \
          = '/netapp/database/pdb/remediated/pdb/%s/pdb%s.ent.gz' \
                % (seed_id[1:3], seed_id[0:4])
    if os.path.exists(original_pdb_filename):
      if ohana:
        cmd = 'ln -s %s %s' % (original_pdb_filename, pdb_filename)
        os.system(cmd)
      else:
        cmd = 'cp %s %s.gz' % (original_pdb_filename, pdb_filename)
        os.system(cmd)
        cmd = 'gunzip %s.gz' % pdb_filename
        os.system(cmd)
    else:
      cmd = "wget %s/%s" % (_PDB_website, pdb_filename)
      os.system(cmd)
  pdb_filename_chain_only = os.path.join(single_seqs_dir(), target_id, 
                                          '%s.pdb' % seed_id[0:5])
  if not os.path.exists(pdb_filename_chain_only):
    if seed_id[4:5] == '':
      cmd = "ln -s %s %s" % (pdb_filename, pdb_filename_chain_only)
      os.system(cmd)
    else:
      in_f = open(pdb_filename)
      out_f = open(pdb_filename_chain_only, "w")
      for line in in_f.readlines():
        if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
          if line[21] == seed_id[4:5]:
            out_f.write(line)
        else:
          out_f.write(line)
      out_f.close()
      in_f.close()
        
    
  status, output = commands.getstatusoutput('sequence_from_pdb.py %s %s' 
                                                % (pdb_filename_chain_only,
                                                  seed_id[4:5]))
  true_seq = output.split('\n')[1]
  # Initialize data structures for true structure and model
  mdl = model(env, file=pdb_filename_chain_only)
  mdl2 = model(env, file=model_file)

  # Create a Modeller .ali file for the target sequence aligned with itself
  # and read this into aln.
  aln = create_target_alignment_with_self_for_modeller(alignment_path,
     seed_id, model_dir, target_id, env, true_seq)

  # Superpose and compute GDT_TS
  (superposition,  num_residues, num_within_1A, num_within_2A, num_within_4A,
    num_within_8A, gdt_ts) = \
	 superpose_and_compute_gdt_ts(mdl, mdl2, aln)

  # Write GDT_TS and its components to a CSV file
  result_file = os.path.join(model_dir, "%s_gdt_ts.csv" % seed_id)
  f = open(result_file, "w")
  f.write("NumResidues,NumWithin1,NumWithin2,NumWithin4,NumWithin8,GDT_TS\n")
  f.write("%d,%d,%d,%d,%d,%g\n" % (num_residues, num_within_1A, num_within_2A,
			 num_within_4A, num_within_8A, gdt_ts))
  f.flush()
  f.close()

  # Record the translation and rotation for later use by make_superposed_pdb.py
  rotation_file = os.path.join(model_dir, "%s_rotation" % seed_id)
  f = open(rotation_file, "w")
  f.write("\n".join([",".join([str(x) for x in row]) 
		      for row in superposition.rotation]))
  f.write("\n")
  f.flush()
  f.close()
  translation_file = os.path.join(model_dir, "%s_translation" % seed_id)
  f = open(translation_file, "w")
  f.write(",".join([str(x) for x in superposition.translation]))
  f.write("\n")
  f.flush()
  f.close()

if __name__ == "__main__":
  main()
