#!/usr/bin/python
'''
Input:
1) For each of two seed sequences, seed_X and seed_Y, presumed
  to be alignable:
     seed sequence, ghmm (HMMs for seed family), and shmms
     (HMMs for seed subfamilies)
2) Reference ("true") alignment between the seeds (this should be optional)
    derived from structural alignment between corresponding structures

Output:
1)  seq alignment, each in its own directory, for each
    combination of (seed, ghmm, ghmm-consensus, hmm, hmm-consensus)
    for seed_X vs. same for seed_Y (stored in "work" directory tree
    contained in seed_X_seed_Y_align directory)
2)  seed-seed alignment, also each in its own directory, derived
    from each of the alignments above (stored in seed_X_seed_Y_align
    directory tree)
3) Cline shift score for each seed-seed alignment vs. reference, stored in
   its own file in same dir as alignment (should be optional)


Method:
  create seq-hmm alignments using align2model
  create hmm-hmm alignments using Lobster
  extract seed-seed alignments from the above
'''

from optparse import OptionParser
import os, glob, re, sys, commands, BPG_common.fasta, string
from array import array
from shmm_shmm_lib import *

global verbose
global hmm_dir
global percent_id


# ------------------------------------------------------------------------------

def shmm_file_of_seed(seed_id, shmm_name):
  return os.path.join(hmm_dir(seed_id, percent_id), "%s.mod" % shmm_name)

# ------------------------------------------------------------------------------

def all_shmm_paths(seed_id):
  # We absolutely require a . before the subfamily identifier
  # Otherwise the ghmm, <seed_id>.mod, will be confused with the shmms
  # Since we have various methods of generating shmms from a fixed tree,
  # and since the name of the ghmm can be confused with the name of the seed
  # sequence, we need to keep these separate.
  return glob.glob(os.path.join(hmm_dir(seed_id, percent_id),
         "%s.*.mod" % seed_id))

# ------------------------------------------------------------------------------

def consensus_seq_file_of_hmm(seed_id, hmm_name):
  return os.path.join(hmm_dir(seed_id, percent_id), "%s.fa" % hmm_name)

# ------------------------------------------------------------------------------

def align_seq_to_hmm(output_dir, seq_id, seq_fasta, hmm_id, hmm_modelfile):
  os.chdir(output_dir)
  cmd = "align2model %s -modelfile %s -db %s &> align2model_%s_%s.out" \
    % ( seq_hmm_file_basename_no_ext(seq_id, hmm_id), hmm_modelfile, seq_fasta,
        seq_id, hmm_id)
  if verbose:
    print "Executing %s in %s" % (cmd, output_dir)
  os.system(cmd)

# ------------------------------------------------------------------------------

def make_prof_of_mod( hmm_name, hmm_path ):
  os.system("sam2psi %s -modelfile %s 2> /dev/null " % (hmm_name, hmm_path))
  cmd = "%s -Ckp2Prof %s.ckp -prof %s.prof" % (lobster_cmd(), hmm_name,
     hmm_name)
  os.system(cmd)

# ------------------------------------------------------------------------------

def get_profile_profile_alignment_path( seedX_id, seedY_id, scoring_function, 
                                        hmmX_name, hmmY_name ):
  combo_id = "%s_%s" % (hmmX_name, hmmY_name)
  return os.path.join(hmmhmm_work_dir(seedX_id, seedY_id, 
                                        hmmX_name, hmmY_name, 0,
                                        scoring_function),
                    "%s_%s_profile_profile.afa" % (combo_id, scoring_function))

# ------------------------------------------------------------------------------

def align_hmm_to_hmm(seedX_id, seedY_id, scoring_function, 
                      hmmX_name, hmmX_path, hmmY_name, hmmY_path):
  profile_profile_alignment_path = \
    get_profile_profile_alignment_path(seedX_id, seedY_id, scoring_function, 
                                        hmmX_name, hmmY_name)
  os.chdir(hmmhmm_work_dir(seedX_id, seedY_id, hmmX_name, hmmY_name, 0,
                            scoring_function))
  make_prof_of_mod(hmmX_name, hmmX_path)
  make_prof_of_mod(hmmY_name, hmmY_path)
  score_output = "%s_%s_%s.out" \
    % (hmmX_name, hmmY_name, scoring_function)
  cmd = "%s -ProfProf -prof1 %s.prof -prof2 %s.prof " \
        % (lobster_cmd(), hmmX_name, hmmY_name) \
        + "-colpair %s -bounds GG -gapstyle simple -out %s > %s" \
        % (scoring_function, profile_profile_alignment_path, score_output)
  print "Aligning profile %s to %s with %s scoring" % (hmmX_name, hmmY_name,
                                              scoring_function)
  os.system(cmd)

# ------------------------------------------------------------------------------

def extract_seed_alignment_from_input_msa(seed_id):
  input_msa = os.path.join(cur_dir(), seed_id, "%s.a2m" % seed_id)
  homolog_seqs = BPG_common.fasta.ReadSequencesList(input_msa)
  for id, seq in homolog_seqs:
    if id[0:len(seed_id)] == seed_id:
      return seq
  print "Error: could not find %s among the %s records in input MSA %s" \
          % (seed_id, len(homolog_seqs), input_msa)
  return ""

# ------------------------------------------------------------------------------

def align_seed_to_own_ghmm(seedX_id, seedY_id, seed_id):
  align_seed_to_ghmm(seedX_id, seedY_id, seed_id, seed_id)

# ------------------------------------------------------------------------------

def align_seed_to_ghmm(seedX_id, seedY_id, seq_seed_id, ghmm_seed_id):
  ghmm_file = ghmm_file_of_seed(ghmm_seed_id)
  seed_fa = os.path.join(single_seqs_dir(), 
                        "%s/%s.fa" % (seq_seed_id, seq_seed_id))
  print "Aligning seed %s to ghmm for %s" % (seq_seed_id, ghmm_seed_id)
  align_seq_to_hmm(hmm_work_dir(seedX_id, seedY_id, ghmm_id(ghmm_seed_id)), 
                  seq_seed_id, seed_fa, ghmm_id(ghmm_seed_id), ghmm_file)

# ------------------------------------------------------------------------------

def align_seed_to_shmm(seedX_id, seedY_id, seed_id, seed_id_of_shmm, shmm_name):
  shmm_file = shmm_file_of_seed(seed_id_of_shmm, shmm_name)
  seed_fa = os.path.join(single_seqs_dir(), "%s/%s.fa" % (seed_id, seed_id))
  align_seq_to_hmm(hmm_work_dir(seedX_id, seedY_id, shmm_name), 
                    seed_id, seed_fa, shmm_name, shmm_file)

# ------------------------------------------------------------------------------

def align_seed_to_shmms(seedX_id, seedY_id, seed_id, seed_id_of_shmms, 
                        shmm_names):
  for shmm_name in shmm_names:
    if verbose:
      print "Aligning seed %s to shmm %s" % (seed_id, shmm_name)
    align_seed_to_shmm(seedX_id, seedY_id, seed_id, seed_id_of_shmms, shmm_name)

# ------------------------------------------------------------------------------

def align_four_way( old_aligned_X_seq, old_aligned_Z_seq, 
                    old_aligned_W_seq, old_aligned_Y_seq ):
  # Input:
  # Z is the consensus sequence of an HMM,
  # and W is the consensus sequence of an HMM (possibly the same one).
  # The Z-HMM and the W-HMM are aligned to each other
  # (possibly trivially, if they are the same HMM).
  # This alignment is manifest in old_aligned_Z_seq and old_aligned_W_seq,
  # which should have the same number of characters (including gaps).
  # Seed sequence X has been aligned to the Z-HMM.
  # This alignment is manifest in old_aligned_X_seq.
  # The number of uppercase characters and dashes in old_aligned_X_seq
  # should be the same as the number of uppercase characters in
  # old_aligned_Z_seq.
  # Seed sequence Y has been aligned to the W-HMM.
  # This alignment is manifest in old_aligned_Y_seq.
  # The number of uppercase characters and dashes in old_aligned_Y_seq
  # should be the same as the number of uppercase characters in
  # old_aligned_W_seq.
  # Output:
  # aligned_X_seq, aligned_Z_seq, aligned_W_seq, aligned_Y_seq
  # a four-way alignment of X, Z, W, and Y
  aligned_X_seq = array('c')
  aligned_Z_seq = array('c')
  aligned_W_seq = array('c')
  aligned_Y_seq = array('c')
  if len(old_aligned_Z_seq) != len(old_aligned_W_seq):
    print "Error in align_four_way"
    print "length %d of aligned Z seq != length %d of aligned W seq" \
      % (len(old_aligned_Z_seq), len(old_aligned_W_seq))
    return None
  else:
    i = 0   # index in X
    j = 0   # index in aligned Z and W
    k = 0   # index in Y
    def go_past_lowercase_X_characters(i):
      while i < len(old_aligned_X_seq) and old_aligned_X_seq[i].islower():
        aligned_X_seq.fromstring(old_aligned_X_seq[i].upper())
        aligned_Z_seq.fromstring("-")
        aligned_W_seq.fromstring("-")
        aligned_Y_seq.fromstring("-")
        i = i + 1
      return i
    def go_past_lowercase_Y_characters(k):
      while k < len(old_aligned_Y_seq) and old_aligned_Y_seq[k].islower():
        aligned_X_seq.fromstring("-")
        aligned_Z_seq.fromstring("-")
        aligned_W_seq.fromstring("-")
        aligned_Y_seq.fromstring(old_aligned_Y_seq[k].upper())
        k = k + 1
      return k
    i = go_past_lowercase_X_characters(i)
    k = go_past_lowercase_Y_characters(k)
    # At the top of this loop, 
    # old_aligned_X_seq[i] is either uppercase or a dash
    # old_aligned_Y_seq[k] is either uppercase or a dash
    for cZ in old_aligned_Z_seq:
      try:
        if cZ.isupper():
          # An uppercase character in the old aligned Z sequence
          # Put the corresponding character from the old alignment of X to Z
          # in the new aligned X sequence, and advance the index
          if i < len(old_aligned_X_seq):
            cX = old_aligned_X_seq[i]
            if cX.isupper() or cX == '-':
              aligned_X_seq.fromstring(old_aligned_X_seq[i])
              i = i + 1
            else:
              print "Error in align_four_way"
              print "Char #%d, %s, in X is neither uppercase nor a dash" \
                % (i, cX)
              return None
          else:
            # If we have run out of characters in the old alignment of X to Z,
            # insert a gap in the new aligned X sequence
            aligned_X_seq.fromstring('-')
        elif cZ == "-":
          # A gap in the old aligned Z sequence
          # Put a new gap in the new aligned X
          aligned_X_seq.fromstring("-")
        else:
          print "Error in align_four_way"
          print "Char #%d, %s, in Z is neither uppercase nor a dash" % (j, cZ)
          return None
        aligned_Z_seq.fromstring(cZ)
        # Z and W are already aligned to each other
        cW = old_aligned_W_seq[j]
        aligned_W_seq.fromstring(cW)
        # Advance the index in the old aligned Z and W sequences
        j = j + 1
        if cW.isupper():
          # An uppercase character in the old aligned W sequence
          # Put the corresponding character from the old alignment of Y to W
          # in the new aligned Y sequence, and advance the index
          # If we have run out of characters in the old alignment of Y to W,
          # insert a gap in the new aligned Y sequence
          if k < len(old_aligned_Y_seq):
            cY = old_aligned_Y_seq[k]
            if cY.isupper() or cY == '-':
              aligned_Y_seq.fromstring(old_aligned_Y_seq[k])
              k = k + 1
            else:
              print "Error in align_four_way"
              print "Char #%d, %s, in Y is neither uppercase nor a dash" \
                % (k, cY)
              return None
          else:
            # If we have run out of characters in the old alignment of Y to W,
            # insert a gap in the new aligned W sequence
            aligned_Y_seq.fromstring('-')
        elif cW == "-":
          # A gap in the old aligned W sequence
          # Put a new gap in the new aligned Y
          aligned_Y_seq.fromstring("-")
        else:
          print "Error in align_four_way"
          print "Char #%d, %s, in W is neither uppercase nor a dash" % (j, cW)
          return None
        # We have now output a single column of the fourway alignment,
        # corresponding to a column of the old Z-W alignment
        # Now we must advance past lowercase characters in the old
        # X-Z and W-Y alignments, 
        # to maintain the invariant at the top of the loop
        i = go_past_lowercase_X_characters(i)
        k = go_past_lowercase_Y_characters(k)
      except IndexError:
        print "IndexError in align_four_way"
        print "i: %d j: %d k: %d" % (i, j, k)
        print "old_aligned_X_seq (%d chars):" % len(old_aligned_X_seq)
        print old_aligned_X_seq
        if i < len(old_aligned_X_seq):
          print "old_aligned_X_seq[%d]: %s" % (i, old_aligned_X_seq[i])
        print "old_aligned_Z_seq (%d chars):" % len(old_aligned_Z_seq)
        print old_aligned_Z_seq
        if j < len(old_aligned_Z_seq):
          print "old_aligned_Z_seq[%d]: %s" % (j, old_aligned_Z_seq[j])
        print "old_aligned_W_seq (%d chars):" % len(old_aligned_W_seq)
        print old_aligned_W_seq
        if j < len(old_aligned_W_seq):
          print "old_aligned_W_seq[%d]: %s" % (j, old_aligned_W_seq[j])
        print "old_aligned_Y_seq (%d chars):" % len(old_aligned_Y_seq)
        print old_aligned_Y_seq
        if k < len(old_aligned_Y_seq):
          print "old_aligned_Y_seq[%d]: %s" % (k, old_aligned_Y_seq[k])
        return None
    # At the end of this loop, we should have used up the characters in 
    # all four sequences

    if i == len(old_aligned_X_seq) and j == len(old_aligned_Z_seq) \
      and j == len(old_aligned_W_seq) and k == len(old_aligned_Y_seq):
      return (aligned_X_seq.tostring(), aligned_Z_seq.tostring(), 
              aligned_W_seq.tostring(), aligned_Y_seq.tostring())
    else:
      print "Error in align_four_way"
      print "Used %d of %d characters of X" % (i, len(old_aligned_X_seq))
      print old_aligned_X_seq
      print "Used %d of %d characters of Z" % (j, len(old_aligned_Z_seq))
      print old_aligned_Z_seq
      print "Used %d of %d characters of W" % (j, len(old_aligned_W_seq))
      print old_aligned_W_seq
      print "Used %d of %d characters of Y" % (k, len(old_aligned_Y_seq))
      print old_aligned_Y_seq
    return None
  
# ------------------------------------------------------------------------------

def align_three_way( old_aligned_X_seq, old_aligned_Z_seq, old_aligned_Y_seq ):
  alignment = align_four_way( old_aligned_X_seq, old_aligned_Z_seq,
                              old_aligned_Z_seq, old_aligned_Y_seq )
  if alignment != None:
    aligned_X_seq, aligned_Z_seq, seq, aligned_Y_seq = alignment
    return (aligned_X_seq, aligned_Z_seq, aligned_Y_seq)
  return None

# ------------------------------------------------------------------------------

def write_fourway(fourway, original_alignment_file, fourway_alignment_file,
                  X_header, Z_header, W_header, Y_header):
  aligned_X_seq, aligned_Z_seq, aligned_W_seq, aligned_Y_seq = fourway
  fourway_fp = open(fourway_alignment_file, "w")
  original_fp = open(original_alignment_file, "w")
  def both_write(str):
    fourway_fp.write(str)
    original_fp.write(str)
  both_write(">%s\n" % X_header)
  both_write("%s\n" % aligned_X_seq)
  fourway_fp.write(">%s\n" % Z_header)
  fourway_fp.write("%s\n" % aligned_Z_seq)
  fourway_fp.write(">%s\n" % W_header)
  fourway_fp.write("%s\n" % aligned_W_seq)
  both_write(">%s\n" % Y_header)
  both_write("%s\n" % aligned_Y_seq)
  fourway_fp.close()
  original_fp.close()
  if verbose:
    print "wrote %s" % original_alignment_file

# ------------------------------------------------------------------------------

def write_threeway(threeway, original_alignment_file, threeway_alignment_file,
                    X_header, Z_header, Y_header):
  aligned_X_seq, aligned_Z_seq, aligned_Y_seq = threeway
  threeway_fp = open(threeway_alignment_file, "w")
  original_fp = open(original_alignment_file, "w")
  def both_write(str):
    threeway_fp.write(str)
    original_fp.write(str)
  both_write(">%s\n" % X_header)
  both_write("%s\n" % aligned_X_seq)
  threeway_fp.write(">%s\n" % Z_header)
  threeway_fp.write("%s\n" % aligned_Z_seq)
  both_write(">%s\n" % Y_header)
  both_write("%s\n" % aligned_Y_seq)
  threeway_fp.close()
  original_fp.close()
  if verbose:
    print "wrote %s" % original_alignment_file

# ------------------------------------------------------------------------------

def compute_sp_cs_score(test_align_file, reference_align_file, out_file):
  os.system("%s -compareMSA %s -ref %s >& %s" \
    % (lobster_cmd(), test_align_file, reference_align_file, out_file))
  
# ------------------------------------------------------------------------------

def read_sp_cs_score(file):
  sp_score_pattern = re.compile("SP=([\-0-9\.eE]*)")
  cs_score_pattern = re.compile("CS=([\-0-9\.eE]*)")
  try:
    line = open(file).read()
    m = re.search(sp_score_pattern, line)
    sp_score = float(m.group(1))
    m = re.search(cs_score_pattern, line)
    cs_score = float(m.group(1))
  except ValueError:
    print "Error in read_sp_cs_score(%s)" % file
    print "line: ", line
    print "Exiting..."
    sys.exit(0)
  return sp_score, cs_score
  
# ------------------------------------------------------------------------------

def write_fourway_and_score( fourway, seedX_id, Z_id, W_id, seedY_id,
                              scoring_function ):
  id = "%s_%s_%s_%s_%s" % (seedX_id, seedY_id, Z_id, W_id, scoring_function )
  shmmX_name = Z_id.split("_")[0]
  shmmY_name = W_id.split("_")[0]
  original_align_file = os.path.join(hmmhmm_align_dir(seedX_id, seedY_id,
                                                      shmmX_name, shmmY_name, 
                                                      0, scoring_function),
                                      "%s_original.afa" % id)
  fourway_align_file = os.path.join(hmmhmm_work_dir(seedX_id, seedY_id, 
                                                    shmmX_name, shmmY_name, 
                                                    0, scoring_function),
                                    "%s_fourway.afa" % id)
  write_fourway(fourway, original_align_file, fourway_align_file,
              seedX_id, Z_id, W_id, seedY_id)
  # compute the CS score
  cs_output_file = os.path.join(hmmhmm_align_dir(seedX_id, seedY_id,
                                                  shmmX_name, shmmY_name, 
                                                  0, scoring_function), 
                                "%s_CS.out" % id)
  compute_sp_cs_score(original_align_file, 
    get_reference_alignment_path("%s_%s" % (seedX_id, seedY_id)),
    cs_output_file)

# ------------------------------------------------------------------------------

def write_threeway_and_score( threeway, seedX_id, Z_id, seedY_id ):
  id = "%s_%s_%s" % (seedX_id, seedY_id, Z_id)
  shmm_name = Z_id.split("_")[0]
  original_align_file = \
    os.path.join(hmm_align_dir(seedX_id, seedY_id, shmm_name), 
                  "%s_original.afa" % id)
  threeway_align_file = \
    os.path.join(hmm_work_dir(seedX_id, seedY_id, shmm_name), 
                  "%s_threeway.afa" % id)
  write_threeway(threeway, original_align_file, threeway_align_file,
                seedX_id, Z_id, seedY_id)
  # compute the CS score
  cs_output_file = os.path.join(hmm_align_dir(seedX_id, seedY_id, shmm_name), 
                                "%s_CS.out" % id)
  compute_sp_cs_score(original_align_file, 
    get_reference_alignment_path("%s_%s" % (seedX_id, seedY_id)),
    cs_output_file)
  
# ------------------------------------------------------------------------------
def align_and_score_originals_from_X_hmm( seedX_id, seedY_id, hmm_name,
                                            input_aligned_X_seq ):
  YZ_file = os.path.join(hmm_work_dir(seedX_id, seedY_id, hmm_name),
                        seq_hmm_file_basename_with_ext(seedY_id, hmm_name))
  old_aligned_Y_seq = BPG_common.fasta.ReadOneSequence(YZ_file)
  Z_file = consensus_seq_file_of_hmm( seedX_id, hmm_name )
  old_aligned_Z_seq = BPG_common.fasta.ReadOneSequence(Z_file)
  threeway = align_three_way( input_aligned_X_seq, old_aligned_Z_seq,
                              old_aligned_Y_seq )
  if threeway == None:
    print "align_three_way failed for hmm %s" % hmm_name
  else:
    write_threeway_and_score(threeway, seedX_id, hmm_name, seedY_id)
# ------------------------------------------------------------------------------
def align_and_score_originals_from_Y_hmm( seedX_id, seedY_id, hmm_name,
                                          input_aligned_Y_seq ):
  XZ_file = os.path.join(hmm_work_dir(seedX_id, seedY_id, hmm_name),
                        seq_hmm_file_basename_with_ext(seedX_id, hmm_name))
  old_aligned_X_seq = BPG_common.fasta.ReadOneSequence(XZ_file)
  Z_file = consensus_seq_file_of_hmm( seedY_id, hmm_name )
  old_aligned_Z_seq = BPG_common.fasta.ReadOneSequence(Z_file)
  threeway = align_three_way( old_aligned_X_seq, old_aligned_Z_seq,
                              input_aligned_Y_seq )
  if threeway == None:
    print "align_three_way failed for hmm %s" % hmm_name
  else:
    write_threeway_and_score(threeway, seedX_id, hmm_name, seedY_id)

# ------------------------------------------------------------------------------

def implied_other_aligned_seq(aligned_seq, other_unaligned_seq):
  if len(aligned_seq) < len(other_unaligned_seq):
    print "Error in implied_other_aligned_seq"
    print "Length of aligned_seq: %d < length of other_unaligned_seq: %d" \
      % (len(aligned_seq), len(other_unaligned_seq))
    return ""
  other_aligned_seq = array('c')
  i = 0
  j = 0
  def go_past_lowercase_characters(i):
    while i < len(aligned_seq) and aligned_seq[i].islower():
      other_aligned_seq.fromstring('-')
      i = i + 1
    return i
  i = go_past_lowercase_characters(i)
  # At the top of this loop, aligned_seq[i] is either uppercase or a dash
  while i < len(aligned_seq):
    try:
      if aligned_seq[i].isupper():
        # An uppercase character in aligned_seq
        # Match it with the corresponding character from other_unaligned_seq
        if other_unaligned_seq[j].isupper():
          other_aligned_seq.fromstring(other_unaligned_seq[j])
          i = i + 1
          j = j + 1
        else:
          print "Error in implied_other_aligned_seq"
          print "Char #%d, %s, of other_unaligned_seq is not uppercase" \
            % (j, other_unaligned_seq[j])
          return ""
      elif aligned_seq[i] == '-':
        # A dash in aligned_seq
        # Insert the corresponding character from other_unaligned_seq
        if other_unaligned_seq[j].isupper():
          other_aligned_seq.fromstring(other_unaligned_seq[j])
          i = i + 1
          j = j + 1
        else:
          print "Error in implied_other_aligned_seq"
          print "Char #%d, %s, of other_unaligned_seq is not uppercase" \
            % (j, other_unaligned_seq[j])
          return ""
      else:
        print "Error in implied_other_aligned_seq"
        print "Char #%d, %s, of aligned_seq is neither uppercase nor a dash" \
          % (i, c)
        return ""
    except IndexError:
      print "IndexError in implied_other_aligned_seq"
      print "i: %d; j: %d" % (i, j)
      print "aligned_seq (%d chars):" % len(aligned_seq)
      print aligned_seq
      if i < len(aligned_seq):
        print "aligned_seq[%d]: %s" % (i, aligned_seq[i])
      print "other_unaligned_seq (%d chars):" % len(other_unaligned_seq)
      print other_unaligned_seq
      if j < len(other_unaligned_seq):
        print "other_unaligned_seq[%d]: %s" % (j, other_unaligned_seq[j])
      return ""
    i = go_past_lowercase_characters(i)
  # At the end of the loop, we should have used up exactly all the characters
  # in aligned_seq and other_unaligned_seq
  if i == len(aligned_seq) and j == len(other_unaligned_seq):
    return other_aligned_seq.tostring()
  else:
    print "Error in implied_other_aligned_seq"
    print "Used %d of %d characters of aligned_seq" % (i, len(aligned_seq))
    print aligned_seq
    print "Used %d of %d characters of other_unaligned_seq" \
      % (j, len(other_unaligned_seq))
    print other_unaligned_seq
  return ""

# ------------------------------------------------------------------------------

def align_and_score_originals_from_Xconsensus_seq_Yhmm( seedX_id, hmmX_name,
                                                        seedY_id, hmmY_name,
                                                        input_aligned_X_seq,
                                                        input_aligned_Y_seq):
  ZW_file = os.path.join(hmmhmm_work_dir(seedX_id, seedY_id, 
                                            hmmX_name, hmmY_name, -1, ""),
                        seq_hmm_file_basename_with_ext(hmmX_name, hmmY_name))
  single_aligned_Z_seq = BPG_common.fasta.ReadOneSequence(ZW_file)
  # The ZW file contains the alignment of a sequence Z,
  # the profile consensus sequence for the HMM named hmmX_name,
  # to an HMM, the one named hmmY_name
  # This implies an alignment with the profile consensus sequence W
  # for the HMM named hmmY_name
  # We make this alignment manifest in old_aligned_W_seq
  W_file = consensus_seq_file_of_hmm( seedY_id, hmmY_name )
  unaligned_W_seq = BPG_common.fasta.ReadOneSequence(W_file)
  old_aligned_W_seq = implied_other_aligned_seq(single_aligned_Z_seq,
                                                unaligned_W_seq)
  # Now dashes were inserted into old_aligned_W_seq
  # corresponding to the lowercase characters of single_aligned_Z_seq
  # so we can uppercase them
  old_aligned_Z_seq = string.upper(single_aligned_Z_seq)
  fourway = align_four_way( input_aligned_X_seq, old_aligned_Z_seq.upper(),
                            old_aligned_W_seq, input_aligned_Y_seq )

  if fourway == None:
    print "align_four_way failed for consensus sequence %s to profile %s" \
      % (hmmX_name, hmmY_name)
  else:
    write_fourway_and_score(fourway, seedX_id, hmmX_name + "_seq", 
                            hmmY_name + "_hmm", seedY_id, "")

# ------------------------------------------------------------------------------

def align_and_score_originals_from_Yconsensus_seq_Xhmm( seedX_id, hmmX_name,
                                                        seedY_id, hmmY_name,
                                                        input_aligned_X_seq,
                                                        input_aligned_Y_seq):
  ZW_file = os.path.join(hmmhmm_work_dir(seedX_id, seedY_id, 
                                            hmmX_name, hmmY_name, 1, ""),
                        seq_hmm_file_basename_with_ext(hmmY_name, hmmX_name))
  single_aligned_W_seq = BPG_common.fasta.ReadOneSequence(ZW_file)
  # The ZW file contains the alignment of a sequence W,
  # the profile consensus sequence for the HMM named hmmY_name,
  # to an HMM, the one named hmmX_name
  # This implies an alignment with the profile consensus sequence Z
  # for the HMM named hmmX_name
  # We make this alignment manifest in old_aligned_Z_seq
  Z_file = consensus_seq_file_of_hmm( seedX_id, hmmX_name )
  unaligned_Z_seq = BPG_common.fasta.ReadOneSequence(Z_file)
  old_aligned_Z_seq = implied_other_aligned_seq(single_aligned_W_seq,
                                                unaligned_Z_seq)
  # Now dashes were inserted into old_aligned_Z_seq
  # corresponding to the lowercase characters of single_aligned_W_seq
  # so we can uppercase them
  old_aligned_W_seq = string.upper(single_aligned_W_seq)
  fourway = align_four_way( input_aligned_X_seq, old_aligned_Z_seq,
                            old_aligned_W_seq.upper(), input_aligned_Y_seq )

  if fourway == None:
    print "align_four_way failed for consensus sequence %s to profile %s" \
      % (hmmY_name, hmmX_name)
  else:
    write_fourway_and_score(fourway, seedX_id, hmmY_name + "_seq", 
                            hmmX_name + "_hmm", seedY_id, "")

def align_and_score_originals_from_hmm_hmm(seedX_id, seedY_id,
                hmmX_name, hmmX_file, hmmY_name, hmmY_file, scoring_function,
                input_aligned_X_seq, input_aligned_Y_seq):
  ZW_file = get_profile_profile_alignment_path(seedX_id, seedY_id,
                            scoring_function, hmmX_name, hmmY_name)
  if (os.path.exists(ZW_file)):
    prof_prof_seqs = BPG_common.fasta.ReadSequencesList(ZW_file)
    old_aligned_Z_seq = ""
    old_aligned_W_seq = ""
    for header, profile_consensus_seq in prof_prof_seqs:
      if len(header) >= len(seedX_id) and header[0:len(seedX_id)] == seedX_id:
        old_aligned_Z_seq = profile_consensus_seq
      elif len(header) >= len(seedY_id) and header[0:len(seedY_id)] == seedY_id:
        old_aligned_W_seq = profile_consensus_seq
      else: 
        print "Error: found extra sequence header %s in %s" % (header, ZW_file)
    fourway = align_four_way( input_aligned_X_seq, old_aligned_Z_seq,
                              old_aligned_W_seq, input_aligned_Y_seq )

    if fourway == None:
      print "align_four_way failed for profile %s to profile %s" \
        % (hmmX_name, hmmY_name)
    else:
      write_fourway_and_score(fourway, seedX_id, hmmX_name, hmmY_name, 
                              seedY_id, scoring_function)
  else:
    print "Error in align_and_score_originals_from_hmm"
    print "File %s does not exist" % ZW_file

# ------------------------------------------------------------------------------


def main():
  # parse command line options
  opt_parser = OptionParser()
  opt_parser.add_option("--pairs", dest="pairs", default="",
                help="comma-separated list of pairs of pdbids seedXid_seedYid")
  opt_parser.add_option("-s", "--use_sciphy_subfams", dest="use_sciphy_subfams",
                        action="store_true", default=False,
               help="Whether to use subfamilies and HMMs generated by SCI-PHY")
  opt_parser.add_option("--use_kerfinfoshare_hmms",
                        dest="use_kerfinfoshare_hmms",
                        action="store_true", default=True,
                 help="Whether to use info-share hmms for kerf subfamilies")
  opt_parser.add_option("--use_kerfw05_hmms",
                        dest="use_kerfinfoshare_hmms",
                        action="store_false", default=True,
                        help="Whether to use w0.5 hmms for kerf subfamilies")
  opt_parser.add_option("-k", "--use_kerf_subfams", dest="use_sciphy_subfams",
                        action="store_false", default=False,
                        help="Whether to use subfamilies defined by %id")
  opt_parser.add_option("-p", "--percentid", type="int", dest="percentid", 
                        default=20,
                        help="Minimum %id at which subfamilies were cut")
  opt_parser.add_option("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Whether to print verbose output")
  opt_parser.add_option("-q", "--quiet", dest="verbose",
                        action="store_false", default=False,
                        help="Whether to suppress verbose output")
  opt_parser.add_option("--profile_profile_scoring_function",
                    dest="scoring_function", default="YL",
                    help="Scoring function for profile-profile alignment")
  (options, args) = opt_parser.parse_args()

  global verbose
  global hmm_dir
  global percent_id
  percent_id = options.percentid
  verbose = options.verbose
  scoring_function = options.scoring_function
  # Split up your list of comma separated pairs into an array of pairs.
  # If you want to use distributed computing, you might just submit one pair
  # at a time, and do it multiple times rather than a big batch list of pairs.
  pairs = string.split(options.pairs, ',')
  if pairs == []:
    os.chdir(pair_dir())
    pairs = glob.glob("*")
  if options.use_sciphy_subfams:
    print "Using SCI-PHY subfamilies and hmms"
    hmm_dir = sciphy_hmm_dir
  else:
    if options.use_kerfinfoshare_hmms:
      hmm_dir = kerf_infoshare_hmm_dir
    else:
      hmm_dir = kerf_w05_hmm_dir
  verbose = options.verbose
  for pair in pairs:
    seedX_id, seedY_id = pair.split("_")
    print "Working on pair %s %s" % (seedX_id, seedY_id)
    print "Reference alignment in %s" % get_reference_alignment_path(pair)

    # -----------------------------------------------------------------------
    # Create directories to hold both pairwise seed-seed alignments
    # (the final alignments between "originals") and the seq-HMM and
    # HMM-HMM alignments needed to generate those seed-seed alignments
    # (the "work" directory tree)
    # -----------------------------------------------------------------------
    make_dir_exist(align_dir(seedX_id, seedY_id))
    print "Will put new pairwise alignments under %s" \
      % align_dir(seedX_id, seedY_id)
    make_dir_exist(work_dir(seedX_id, seedY_id))
    print "Will put other new output files under %s" \
      % work_dir(seedX_id, seedY_id)

    # Make alignment/modelling directories for ghmms.

    # GHMM_X_SEED_Y
    make_dir_exist(hmm_align_dir(seedX_id, seedY_id, ghmm_id(seedX_id)))
    # GHMM_X_GHMM_Y and associated HMM_SEQ combos
    for cons in range(-1,2):
      make_dir_exist(hmmhmm_align_dir(seedX_id, seedY_id,
                            ghmm_id(seedX_id), ghmm_id(seedY_id), cons,
                            scoring_function))
    # GHMM_Y_SEED_X
    make_dir_exist(hmm_align_dir(seedX_id, seedY_id, ghmm_id(seedY_id)))
    # GHMM_Y_GHMM_X  and associated HMM_SEQ combos (SYMLINK)
    make_hmmhmm_symlinks(seedX_id, seedY_id, ghmm_id(seedX_id),
                          ghmm_id(seedY_id), scoring_function, hmm_align_dir)

    # Make parallel work directories for ghmms.

    # GHMM_X_SEED_Y
    make_dir_exist(hmm_work_dir(seedX_id, seedY_id, ghmm_id(seedX_id)))
    # GHMM_X_GHMM_Y and associated HMM_SEQ combos
    for cons in range(-1,2):
      make_dir_exist(hmmhmm_work_dir(seedX_id, seedY_id,
                            ghmm_id(seedX_id), ghmm_id(seedY_id), cons,
                            scoring_function))
    # GHMM_Y_SEED_Y
    make_dir_exist(hmm_work_dir(seedX_id, seedY_id, ghmm_id(seedY_id)))
    # GHMM_Y_GHMM_X  and associated HMM_SEQ combos (SYMLINK)
    make_hmmhmm_symlinks(seedX_id, seedY_id, ghmm_id(seedX_id),
                          ghmm_id(seedY_id), scoring_function, hmm_work_dir)

    # Make alignment/modelling directories for shmms
 
    # locate all the shmms for X
    shmmX_paths = all_shmm_paths(seedX_id)
    print "Found %d shmms for %s" % (len(shmmX_paths), seedX_id)
    shmmX_names = []
    for shmmX_path in shmmX_paths:
      shmmX_names.append(os.path.splitext(os.path.split(shmmX_path)[1])[0])
    # locate all the shmms for Y
    shmmY_paths = all_shmm_paths(seedY_id)
    print "Found %d shmms for %s" % (len(shmmY_paths), seedY_id)
    shmmY_names = []
    for shmmY_path in shmmY_paths:
      shmmY_names.append(os.path.splitext(os.path.split(shmmY_path)[1])[0])


    # make alignment and work directories for each shmm of X
    for shmmX_name in shmmX_names:

      # SHMM_X_SEED_Y
      make_dir_exist(hmm_align_dir(seedX_id, seedY_id, shmmX_name))
      make_dir_exist(hmm_work_dir(seedX_id, seedY_id, shmmX_name))

      # SHMM_X_GHMM_Y  and  GHMM_Y_SHMM_X, and HMM_SEQ combos (SYMLINK)
      for cons in range(-1,2):
        make_dir_exist(hmmhmm_align_dir(seedX_id, seedY_id, shmmX_name,
                                        ghmm_id(seedY_id), cons,
                                        scoring_function))
        make_dir_exist(hmmhmm_work_dir(seedX_id, seedY_id, shmmX_name,
                                        ghmm_id(seedY_id), cons,
                                        scoring_function))
      dest = hmmhmm_align_dir(seedX_id, seedY_id,
                                   ghmm_id(seedY_id), shmmX_name, 0,
                                   scoring_function)
      make_hmmhmm_symlinks(seedX_id, seedY_id, shmmX_name, ghmm_id(seedY_id),
                            scoring_function, hmm_align_dir)
      make_hmmhmm_symlinks(seedX_id, seedY_id, shmmX_name, ghmm_id(seedY_id),
                            scoring_function, hmm_work_dir)

      # each contains an alignment to each shmm of Y
      for shmmY_name in shmmY_names:
        # SHMM_X_SHMM_Y and HMM_SEQ combos
        for cons in range(-1,2):
          make_dir_exist(hmmhmm_align_dir(seedX_id, seedY_id, 
                                          shmmX_name, shmmY_name, cons,
                                          scoring_function))
          make_dir_exist(hmmhmm_work_dir(seedX_id, seedY_id, 
                                          shmmX_name, shmmY_name, cons,
                                          scoring_function))

    # make alignment and work directories for each shmm of Y
    for shmmY_name in shmmY_names:
      # SHMM_Y_SEED_X
      make_dir_exist(hmm_align_dir(seedX_id, seedY_id, shmmY_name))
      make_dir_exist(hmm_work_dir(seedX_id, seedY_id, shmmY_name))
      # GHMM_X_SHMM_Y  and  SHMM_Y_GHMM_X, and HMM_SEQ combos (SYMLINK)
      for cons in range(-1,2):
        make_dir_exist(hmmhmm_align_dir(seedX_id, seedY_id, ghmm_id(seedX_id),
                          shmmY_name, cons, scoring_function))
        make_dir_exist(hmmhmm_work_dir(seedX_id, seedY_id, ghmm_id(seedX_id),
                          shmmY_name, cons, scoring_function))
      make_hmmhmm_symlinks(seedX_id, seedY_id, ghmm_id(seedX_id), shmmY_name,
                            scoring_function, hmm_align_dir)
      make_hmmhmm_symlinks(seedX_id, seedY_id, ghmm_id(seedX_id), shmmY_name,
                            scoring_function, hmm_work_dir)

      # each contains an alignment to each shmm of X
      for shmmX_name in shmmX_names:
        # SHMM_Y_SHMM_X  and HMM_SEQ combos (SYMLINK)
        make_hmmhmm_symlinks(seedX_id, seedY_id, shmmX_name, shmmY_name,
                              scoring_function, hmm_align_dir)
        make_hmmhmm_symlinks(seedX_id, seedY_id, shmmX_name, shmmY_name,
                              scoring_function, hmm_work_dir)

   # ------------------------------------------------------------------
   # Create seed-HMM and HMM-HMM alignments and store in "work" dir
   # ------------------------------------------------------------------

    input_aligned_X_seq = ""
    input_aligned_Y_seq = ""

    print "Extracting aligned seed sequences from input MSA"
    input_aligned_X_seq = extract_seed_alignment_from_input_msa(seedX_id)
    input_aligned_Y_seq = extract_seed_alignment_from_input_msa(seedY_id)

    print "Aligning seeds to the other ghmms and shmms"
    align_seed_to_ghmm(seedX_id, seedY_id, seedX_id, seedY_id)
    align_seed_to_ghmm(seedX_id, seedY_id, seedY_id, seedX_id)
    align_seed_to_shmms(seedX_id, seedY_id, seedX_id, seedY_id, shmmY_names)
    align_seed_to_shmms(seedX_id, seedY_id, seedY_id, seedX_id, shmmX_names)

    # Profile-profile alignments
    # GHMM_X_GHMM_Y
    print "Aligning %s to %s using %s profile-profile alignment" \
      % (ghmm_id(seedX_id), ghmm_id(seedY_id), scoring_function)
    align_hmm_to_hmm( seedX_id, seedY_id, scoring_function, 
                      ghmm_id(seedX_id), ghmm_file_of_seed(seedX_id),
                      ghmm_id(seedY_id), ghmm_file_of_seed(seedY_id))
    # SHMM_X_GHMM_Y
    print "Aligning shmms to the other ghmm by %s profile-profile alignment" \
      % scoring_function
    for shmmX_name in shmmX_names:
      if verbose:
        print "Aligning %s to %s" % (shmmX_name, ghmm_id(seedY_id))
      align_hmm_to_hmm( seedX_id, seedY_id, scoring_function,
                      shmmX_name, shmm_file_of_seed(seedX_id, shmmX_name),
                      ghmm_id(seedY_id), ghmm_file_of_seed(seedY_id))
    # GHMM_X_SHMM_Y
    for shmmY_name in shmmY_names:
      if verbose:
        print "Aligning %s to %s" % (ghmm_id(seedX_id), shmmY_name)
      align_hmm_to_hmm( seedX_id, seedY_id, scoring_function,
                      ghmm_id(seedX_id), ghmm_file_of_seed(seedX_id),
                      shmmY_name, shmm_file_of_seed(seedY_id, shmmY_name))
    # SHMM_SHMM
    # Align each of the SHMMs for X to each of the SHMMs for Y
    print "Aligning shmms to shmms by %s profile-profile alignment" \
      % scoring_function
    for shmmX_name in shmmX_names:
      for shmmY_name in shmmY_names:
        if verbose:
          print "Aligning %s to %s" % (shmmX_name, shmmY_name)
        align_hmm_to_hmm( seedX_id, seedY_id,
                          scoring_function,
                          shmmX_name, shmm_file_of_seed(seedX_id, shmmX_name),
                          shmmY_name, shmm_file_of_seed(seedY_id, shmmY_name))
                

    # Get the profile consensus sequence of each ghmm
    # from the ghmm-ghmm profile-profile alignment
    # Align the profile consensus sequence of each ghmm to the other ghmm
    print "Getting the profile consensus sequences of the ghmms"
    # the following 6 lines used by lots of the script. PUT EARLIER
    #  so can easily "comment in" when commenting out chunks
    ghmm_ghmm_file = \
      get_profile_profile_alignment_path(seedX_id, seedY_id, 
                scoring_function, ghmm_id(seedX_id), ghmm_id(seedY_id))
    ghmmghmm_seqs = BPG_common.fasta.ReadSequencesList(ghmm_ghmm_file)
    ghmmX_fa = consensus_seq_file_of_hmm(seedX_id, ghmm_id(seedX_id))
    ghmmY_fa = consensus_seq_file_of_hmm(seedY_id, ghmm_id(seedY_id))
    for header, seq in ghmmghmm_seqs:
      if len(header) >= len(seedX_id) and header[0:len(seedX_id)] == seedX_id:
        fp = open(ghmmX_fa, "w")
        fp.write(">%s\n" % ghmm_id(seedX_id))
        fp.write("%s\n" % string.upper(string.replace(seq, "-", "")))
        fp.close()
      elif len(header) >= len(seedY_id) and header[0:len(seedY_id)] == seedY_id:
        fp = open(ghmmY_fa, "w")
        fp.write(">%s\n" % ghmm_id(seedY_id))
        fp.write("%s\n" % string.upper(string.replace(seq, "-", "")))
        fp.close()
      else:
        print "Error: extra header %s in %s" % (header, ghmm_ghmm_file)
    print "Aligning profile consensus sequence of each ghmm to the other"
    # GHMM_SEQ_X_GHMM_Y
    align_seq_to_hmm( hmmhmm_work_dir(seedX_id, seedY_id, 
                            ghmm_id(seedX_id), ghmm_id(seedY_id), -1, ""), 
                      ghmm_id(seedX_id), ghmmX_fa,
                      ghmm_id(seedY_id), ghmm_file_of_seed(seedY_id))
    # GHMM_X_GHMM_SEQ_Y
    align_seq_to_hmm( hmmhmm_work_dir(seedX_id, seedY_id,
                            ghmm_id(seedX_id), ghmm_id(seedY_id), 1, ""), 
                      ghmm_id(seedY_id), ghmmY_fa,
                      ghmm_id(seedX_id), ghmm_file_of_seed(seedX_id))

    # Get the profile consensus sequence for each shmm
    # from one of the profile-profile alignments
    # Align it to the other ghmm
    # Also align the profile consensus sequence of the other ghmm to the shmm
    print "Getting profile consensus sequences for each shmm"
    for shmmX_name in shmmX_names:
      if verbose:
        print "Getting profile consensus sequence for %s" % shmmX_name
      ghmmY_shmmX_file = \
        get_profile_profile_alignment_path(seedX_id, seedY_id,
                        scoring_function, shmmX_name, ghmm_id(seedY_id))
      ghmmYshmmX_seqs = BPG_common.fasta.ReadSequencesList(ghmmY_shmmX_file)
      shmmX_fa = consensus_seq_file_of_hmm(seedX_id, shmmX_name)
      for header, seq in ghmmYshmmX_seqs:
        if len(header) < len(seedY_id) or \
            header[0:len(seedY_id)] != seedY_id:
          fp = open(shmmX_fa, "w")
          fp.write(">%s\n" % shmmX_name)
          fp.write("%s\n" % string.upper(string.replace(seq, "-", "")))
          fp.close()
      if verbose:
        print "Aligning it to the ghmm for %s" % seedY_id
      # SHMM_SEQ_X_GHMM_Y
      align_seq_to_hmm(hmmhmm_work_dir(seedX_id, seedY_id, 
                                        shmmX_name, ghmm_id(seedY_id), -1, ""),  
        shmmX_name, shmmX_fa, ghmm_id(seedY_id), ghmm_file_of_seed(seedY_id))
      if verbose:
        print "Aligning profile consensus sequence of ghmm for %s to %s" \
              % (seedY_id, shmmX_name)
      # SHMM_X_GHMM_SEQ_Y
      align_seq_to_hmm(hmmhmm_work_dir(seedX_id, seedY_id, 
                                        shmmX_name, ghmm_id(seedY_id), 1, ""), 
        ghmm_id(seedY_id), ghmmY_fa, shmmX_name, 
        shmm_file_of_seed(seedX_id, shmmX_name))
    for shmmY_name in shmmY_names:
      if verbose:
        print "Getting profile consensus sequence for %s" % shmmY_name
      ghmmX_shmmY_file = \
        get_profile_profile_alignment_path(seedX_id, seedY_id,
                        scoring_function, ghmm_id(seedX_id), shmmY_name)
      ghmmXshmmY_seqs = BPG_common.fasta.ReadSequencesList(ghmmX_shmmY_file)
      shmmY_fa = consensus_seq_file_of_hmm(seedY_id, shmmY_name)
      for header, seq in ghmmXshmmY_seqs:
        if len(header) < len(seedX_id) or \
            header[0:len(seedX_id)] != seedX_id:
          fp = open(shmmY_fa, "w")
          fp.write(">%s\n" % shmmY_name)
          fp.write("%s\n" % string.upper(string.replace(seq, "-", "")))
          fp.close()
      if verbose:
        print "Aligning it to the ghmm for %s" % seedX_id
      # GHMM_X_SHMM_SEQ_Y
      next_dir = hmmhmm_work_dir(seedX_id, seedY_id, 
                                        ghmm_id(seedX_id), shmmY_name, 1, "")
      align_seq_to_hmm(next_dir,
        shmmY_name, shmmY_fa, ghmm_id(seedX_id), ghmm_file_of_seed(seedX_id))
      if verbose:
        print "Aligning profile consensus sequence of ghmm for %s to %s" \
              % (seedX_id, shmmY_name)
      # GHMM_SEQ_X_SHMM_Y
      align_seq_to_hmm(hmmhmm_work_dir(seedX_id, seedY_id, 
                                        ghmm_id(seedX_id), shmmY_name, -1, ""), 
        ghmm_id(seedX_id), ghmmX_fa, shmmY_name, 
        shmm_file_of_seed(seedY_id, shmmY_name))
    # Align the profile consequence of each shmm to the other shmm
    for shmmX_name in shmmX_names:
      for shmmY_name in shmmY_names:
        if verbose:
          print "Aligning profile consensus sequence of %s to %s" \
            % (shmmX_name, shmmY_name)
        shmmX_fa = consensus_seq_file_of_hmm(seedX_id, shmmX_name)
        # SHMM_SEQ_X_SHMM_Y
        align_seq_to_hmm(hmmhmm_work_dir(seedX_id, seedY_id,
                                          shmmX_name, shmmY_name, -1, ""), 
           shmmX_name, shmmX_fa, shmmY_name, 
           shmm_file_of_seed(seedY_id, shmmY_name))
        if verbose:
          print "Aligning profile consensus sequence of %s to %s" \
            % (shmmY_name, shmmX_name)
        shmmY_fa = consensus_seq_file_of_hmm(seedY_id, shmmY_name)
        # SHMM_X_SHMM_SEQ_Y
        align_seq_to_hmm(hmmhmm_work_dir(seedX_id, seedY_id, 
                                            shmmX_name, shmmY_name, 1, ""), 
            shmmY_name, shmmY_fa, shmmX_name, 
            shmm_file_of_seed(seedX_id, shmmX_name))

   # ------------------------------------------------------------------
   # Extract pairwise alignments of seeds ("originals") from alignments
   # of HMMs, and compute CS score for each. Each alignment goes into
   # a separate file in a separate directory. CS score goes into its
   # own file in the same dir as its alignment.
   # ------------------------------------------------------------------

    print "Extracting pairwise alignments and scoring"
    print "Using alignments of a seed sequence to an HMM"
    print "Using alignment of %s to the ghmm for %s" % (seedX_id, seedY_id)
    align_and_score_originals_from_Y_hmm(seedX_id, seedY_id, ghmm_id(seedY_id),
                                          input_aligned_Y_seq)
    print "Using alignment of %s to the ghmm for %s" % (seedY_id, seedX_id)
    align_and_score_originals_from_X_hmm(seedX_id, seedY_id, ghmm_id(seedX_id),
                                          input_aligned_X_seq)
    print "Using alignments of %s to the shmms for %s" % (seedY_id, seedX_id)
    for shmmX_name in shmmX_names:
      if verbose:
        print "Using alignment of %s to %s" % (seedY_id, shmmX_name)
      align_and_score_originals_from_X_hmm(seedX_id, seedY_id, shmmX_name,
                                            input_aligned_X_seq)
    print "Using alignments of %s to the shmms for %s" % (seedX_id, seedY_id)
    for shmmY_name in shmmY_names:
      if verbose:
        print "Using alignment of %s to %s" % (seedX_id, shmmY_name)
      align_and_score_originals_from_Y_hmm(seedX_id, seedY_id, shmmY_name,
                                            input_aligned_Y_seq)

    print "Using alignments of a profile consensus sequence to an HMM"
    print "Consensus sequence of GHMM of %s to GHMM of %s" \
          % (seedX_id, seedY_id)
    align_and_score_originals_from_Xconsensus_seq_Yhmm(seedX_id,
                                ghmm_id(seedX_id), seedY_id, ghmm_id(seedY_id),
                                input_aligned_X_seq, input_aligned_Y_seq)
    print "Consensus sequence of GHMM of %s to GHMM of %s" \
          % (seedY_id, seedX_id)
    align_and_score_originals_from_Yconsensus_seq_Xhmm(seedX_id,
                                ghmm_id(seedX_id), seedY_id, ghmm_id(seedY_id),
                                input_aligned_X_seq, input_aligned_Y_seq)
    print "Consensus sequences of SHMMs of %s to GHMM of %s and " \
          % (seedX_id, seedY_id) + \
          "consensus sequence of GHMM of %s to SHMMs of %s" \
          % (seedY_id, seedX_id)
    for shmmX_name in shmmX_names:
      if verbose:
        print "Consensus sequence of %s to ghmm of %s" % (shmmX_name, seedY_id)
      align_and_score_originals_from_Xconsensus_seq_Yhmm( seedX_id, shmmX_name,
                                                  seedY_id, ghmm_id(seedY_id),
                                                  input_aligned_X_seq,
                                                  input_aligned_Y_seq)
      if verbose:
        print "Consensus sequence of ghmm of %s to %s" % (seedY_id, shmmX_name)
      align_and_score_originals_from_Yconsensus_seq_Xhmm( seedX_id, shmmX_name,
                                                  seedY_id, ghmm_id(seedY_id),
                                                  input_aligned_X_seq,
                                                  input_aligned_Y_seq)
    print "Consensus sequences of SHMMs of %s to GHMM of %s and " \
          % (seedY_id, seedX_id) + \
          "consensus sequence of GHMM of %s to SHMMs of %s" \
          % (seedX_id, seedY_id)
    for shmmY_name in shmmY_names:
      if verbose:
        print "Consensus sequence of %s to ghmm of %s" % (shmmY_name, seedX_id)
      align_and_score_originals_from_Yconsensus_seq_Xhmm( seedX_id,
                                        ghmm_id(seedX_id), seedY_id, shmmY_name,
                                        input_aligned_X_seq,
                                        input_aligned_Y_seq)
      if verbose:
        print "Consensus sequence of ghmm of %s to %s" % (seedX_id, shmmY_name)
      align_and_score_originals_from_Xconsensus_seq_Yhmm( seedX_id,
                                        ghmm_id(seedX_id), seedY_id, shmmY_name,
                                        input_aligned_X_seq,
                                        input_aligned_Y_seq)
    for shmmX_name in shmmX_names:
      for shmmY_name in shmmY_names:
        if verbose:
          print "Consensus sequence of %s to %s" % (shmmX_name, shmmY_name)
        align_and_score_originals_from_Xconsensus_seq_Yhmm( seedX_id,
                                        shmmX_name, seedY_id, shmmY_name,
                                        input_aligned_X_seq,
                                        input_aligned_Y_seq)
        if verbose:
          print "Consensus sequence of %s to %s" % (shmmY_name, shmmX_name)
        align_and_score_originals_from_Yconsensus_seq_Xhmm( seedX_id,
                                        shmmX_name, seedY_id, shmmY_name,
                                        input_aligned_X_seq,
                                        input_aligned_Y_seq)


    print "Using %s profile-profile alignments of two HMMs" \
      % scoring_function
    print "Using alignment of %s to %s" % (ghmm_id(seedX_id), ghmm_id(seedY_id))
    align_and_score_originals_from_hmm_hmm(seedX_id, seedY_id,
                    ghmm_id(seedX_id), ghmm_file_of_seed(seedX_id),
                    ghmm_id(seedY_id), ghmm_file_of_seed(seedY_id),
                    scoring_function, input_aligned_X_seq, input_aligned_Y_seq)
    print "Using alignments of %s to shmms for %s" \
      % (ghmm_id(seedY_id), seedX_id)
    for shmmX_name in shmmX_names:
      align_and_score_originals_from_hmm_hmm(seedX_id, seedY_id,
                    shmmX_name, shmm_file_of_seed(seedX_id, shmmX_name),
                    ghmm_id(seedY_id), ghmm_file_of_seed(seedY_id),
                    scoring_function, input_aligned_X_seq, input_aligned_Y_seq)
    print "Using alignments of %s to shmms for %s" \
      % (ghmm_id(seedX_id), seedY_id)
    for shmmY_name in shmmY_names:
      align_and_score_originals_from_hmm_hmm(seedX_id, seedY_id,
                    ghmm_id(seedX_id), ghmm_file_of_seed(seedX_id),
                    shmmY_name, shmm_file_of_seed(seedY_id, shmmY_name),
                    scoring_function, input_aligned_X_seq, input_aligned_Y_seq)
    print "Using alignments of shmms for %s to shmms for %s" \
      % (seedX_id, seedY_id)
    for shmmX_name in shmmX_names:
      for shmmY_name in shmmY_names:
        align_and_score_originals_from_hmm_hmm(seedX_id, seedY_id,
                    shmmX_name, shmm_file_of_seed(seedX_id, shmmX_name),
                    shmmY_name, shmm_file_of_seed(seedY_id, shmmY_name),
                    scoring_function, input_aligned_X_seq, input_aligned_Y_seq)
    print "align_originals_via_shmms_and_score is done with %s" % pair

if __name__ == "__main__":
  main()
