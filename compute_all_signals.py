#!/usr/bin/env python

import os, sys, re, glob, math, cPickle, tempfile
from matchmaker.shmm_shmm_lib import *
from matchmaker.alignments_lib import *
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62

class SecondaryStructurePrediction:
  def __init__(self):
    self.ss = {}
    self.conf = {}
    self.prob = {}

  def read(self, seed_id):
    f = open(os.path.join(single_seqs_dir(), seed_id, 'vanilla_psipred',
                          "%s.ss2" % seed_id))
    lines = [line.rstrip() for line in f.readlines()]
    f.close()
    for line in lines:
      if len(line) > 0 and line[0] != '#':
        fields = line.split()
        i = int(fields[0])-1
        self.ss[i] = fields[2]
        self.prob[i] = [float(x) for x in fields[3:]]
    f = open(os.path.join(single_seqs_dir(), seed_id, 'vanilla_psipred',
                          "%s.horiz" % seed_id))
    lines = [line.rstrip() for line in f.readlines()]
    f.close()
    i = 0
    for line in lines:
      if line.startswith('Conf:'):
        line = line[5:].strip()
        for c in line:
          self.conf[i] = int(c)
          i += 1

def workdir_of_profile_profile_alignment(alignment_file):
  dir_components = alignment_file.split('/')
  workdir = '/'.join(dir_components[:-4]) + '/work/' + \
            '/'.join(dir_components[-4:-1])
  return workdir

def jensen_shannon_divergence(p, q):
  if len(p) != len(q):
    print "Jensen-Shannon divergence requires same sample space"
    return
  sum = 0.0
  m = {}
  for k in xrange(len(p)):
    m[k] = (p[k] + q[k])/2
    if m[k] > 0.0:
      try:
        sum += (p[k] * math.log( p[k] / m[k] ) \
                + q[k] * math.log( q[k] / m[k] ))/2
      except OverflowError:
        pass
  return sum
    
def hmmsum_score(aa0, aa1, gammas0, gammas1, matrices):
  score = 0.0
  for k in xrange(281):
    score += gammas0[k] * matrices[k][(aa0, aa1)] * gammas1[k]
  return score

def main():
  if len(sys.argv) < 2:
    print "Usage: %s <pair>" % sys.argv[0]
    sys.exit(0)

  pair = sys.argv[1]
  seedX, seedY = pair.split('_')
  YL_pat = re.compile('_YL_')
  hhalignsamss_pat = re.compile('_hhalignsamss_')
  hhalign_score_line_header_pat = re.compile(' No Hit')
  f = open(os.path.join(single_seqs_dir(), seedX, "%s.fa" % seedX))
  seedX_record = SeqIO.parse(f, "fasta").next()
  f.close()
  f = open(os.path.join(single_seqs_dir(), seedY, "%s.fa" % seedY))
  seedY_record = SeqIO.parse(f, "fasta").next()
  f.close()
  seqX = seedX_record.seq.tostring()
  seqY = seedY_record.seq.tostring()
  x_gammas = {}
  y_gammas = {}
  f = open(os.path.join(single_seqs_dir(), seedX, "%s.gca" % seedX))
  lines = f.readlines()
  f.close()
  i = 0
  for line in lines:
    if line.startswith('RESIDUE'):
      i = int(line.split()[1]) - 1
    elif line.startswith('GAMMA'):
      x_gammas[i] = [float(x) for x in line[6:].strip().split()]
  f = open(os.path.join(single_seqs_dir(), seedY, "%s.gca" % seedY))
  lines = f.readlines()
  f.close()
  j = 0
  for line in lines:
    if line.startswith('RESIDUE'):
      j = int(line.split()[1]) - 1
    elif line.startswith('GAMMA'):
      y_gammas[j] = [float(x) for x in line[6:].strip().split()]
  f = open("/clusterfs/ohana/software/HMMSUM/matrices/observed_NS.pkl")
  hmmsum_matrices = cPickle.load(f)
  f.close()
  ss_X = SecondaryStructurePrediction()
  ss_X.read(seedX)
  ss_Y = SecondaryStructurePrediction()
  ss_Y.read(seedY)
  f = open(os.path.join(single_seqs_dir(), seedX, 'sable_output_graph'))
  lines = [line.strip() for line in f.readlines()]
  f.close()
  sa_X = [int(x) for x in lines[4].split()]
  sa_conf_X = [int(lines[5][k]) for k in range(len(lines[5]))]
  f = open(os.path.join(single_seqs_dir(), seedY, 'sable_output_graph'))
  lines = [line.strip() for line in f.readlines()]
  f.close()
  sa_Y = [int(x) for x in lines[4].split()]
  sa_conf_Y = [int(lines[5][k]) for k in range(len(lines[5]))]
  f = open(os.path.join(single_seqs_dir(), seedX, 'rs3_%s.out' % seedX))
  lines = [line.rstrip() for line in f.readlines()]
  f.close()
  phi_X = {}
  psi_X = {}
  for line in lines:
    if len(line) > 0 and line[0] != '>' and line[0] != '#':
      index, aa, asa, phi, psi = line.strip().split()
      i = int(index) - 1
      phi_X[i] = float(phi)
      psi_X[i] = float(psi)
  f = open(os.path.join(single_seqs_dir(), seedY, 'rs3_%s.out' % seedY))
  lines = [line.rstrip() for line in f.readlines()]
  f.close()
  phi_Y = {}
  psi_Y = {}
  for line in lines:
    if len(line) > 0 and line[0] != '>' and line[0] != '#':
      index, aa, asa, phi, psi = line.strip().split()
      i = int(index) - 1
      phi_Y[i] = float(phi)
      psi_Y[i] = float(psi)
  qscore_directory = tempfile.mkdtemp()
  temp_afa_name = os.path.join(qscore_directory, "temp.afa")
  ref_afa_path = get_reference_alignment_path(pair)
  max_q_combined = 0.0
  max_q_developer = 0.0
  max_q_modeler = 0.0
  max_cline_shift = 0.0
  max_total_column = 0.0
  max_num_aligned_pairs = 0
  max_yl_score = -10000.0
  max_hhalign_score = -10000.0
  max_hhalign_ss = -10000.0
  max_hhalign_cols = 0
  max_percent_id = 0.0
  max_mean_blosum62_score = -10000.0
  max_mean_hmmsum_score = -10000.0
  max_fraction_ss_agreement = 0.0
  max_weighted_fraction_ss_agreement = 0.0
  max_total_secondary_structure_weights = 0.0
  max_mean_ss_dot_product = 0.0
  min_mean_solvent_accessibility_difference = 100
  min_weighted_mean_solvent_accessibility_difference = 10000
  max_total_solvent_accessibility_weights = 0
  min_mean_phi_difference = 180.0
  min_mean_psi_difference = 180.0
  sys.stdout.write("Alignment,Qcombined,")
  sys.stdout.write("Qdeveloper,Qmodeler,ClineShift,")
  sys.stdout.write("NumAlignedPairs,")
  sys.stdout.write("YL,HHalignScore,HHalignSS,HHalignCols,")
  sys.stdout.write("PercentId,MeanBLOSUM62Score,")
  sys.stdout.write("MeanHMMSUMScore,")
  sys.stdout.write("FractionPredictedSecondaryStructureAgreement,")
  sys.stdout.write(
              "ConfidenceWeightedFractionPredictedSecondaryStructureAgreement,")
  sys.stdout.write("TotalConfidenceWeightsPredictedSecondaryStructure,")
  sys.stdout.write("MeanSecondaryStructureDotProduct,")
  sys.stdout.write("MeanSolventAccessibilityDifference,")
  sys.stdout.write("ConfidenceWeightedMeanSolventAccessibilityDifference,")
  sys.stdout.write("TotalConfidenceWeightsPredictedSolventAccessibility,")
  sys.stdout.write("MeanPhiDifference,")
  # Last feature - no comma
  sys.stdout.write("MeanPsiDifference")
  sys.stdout.write("\n")
  qf = open(os.path.join(cur_dir(), "results", "%s_align" % pair,
            "%s_q_combined_scores.csv" % pair))
  for qf_line in qf.readlines():
    alignment_name, q_combined = qf_line.strip().split(',')
    q_combined = float(q_combined)
    if q_combined > max_q_combined:
      max_q_combined = q_combined
    alignment_file = os.path.join(alignment_dir(alignment_name),
                                  alignment_name)
    sys.stdout.write("%s," % qf_line.rstrip())
    ((header1, seq1), (header2, seq2)) = read_alignment_file(alignment_file)
    f = open(temp_afa_name, "w")
    f.write(">%s\n" % header1)
    f.write("%s\n" % seq1)
    f.write(">%s\n" % header2)
    f.write("%s\n" % seq2)
    f.close()
    cmd = "qscore -test %s -ref %s" % (temp_afa_name, ref_afa_path)
    status, output = commands.getstatusoutput(cmd)
    scores = output.strip().split(';')
    q_developer = float(scores[2].split('=')[1])
    q_modeler = float(scores[3].split('=')[1])
    cline_shift = float(scores[4].split('=')[1])
    sys.stdout.write("%g,%g,%g," % (q_developer, q_modeler, cline_shift))
    if q_developer > max_q_developer:
      max_q_developer = q_developer
    if q_modeler > max_q_modeler:
      max_q_modeler = q_modeler
    if cline_shift > max_cline_shift:
      max_cline_shift = cline_shift
    aln = Alignment()
    aln.readFromFile(alignment_file)
    num_aligned_pairs = len(aln.pairsNoGaps)
    if num_aligned_pairs > max_num_aligned_pairs:
      max_num_aligned_pairs = num_aligned_pairs
    sys.stdout.write("%d," % num_aligned_pairs)
    if YL_pat.search(alignment_file):
      workdir = workdir_of_profile_profile_alignment(alignment_file)
      f = open(glob.glob(os.path.join(workdir, '*_YL.out'))[0])
      YL_score = f.read().strip().split()[2]
      f.close()
      sys.stdout.write("%s,,,," % YL_score)
      YL_score = float(YL_score)
      if YL_score > max_yl_score:
        max_yl_score = YL_score
    elif hhalignsamss_pat.search(alignment_file):
      workdir = workdir_of_profile_profile_alignment(alignment_file)
      f = open(glob.glob(os.path.join(workdir, '*.align'))[0])
      lines = [line.rstrip() for line in f.readlines()]
      f.close()
      found_score_line_header = False
      for line in lines:
        if line.startswith(' No Hit'):
          found_score_line_header = True
        elif found_score_line_header:
          hit_number, hit, prob, evalue, pvalue, score, ss, cols, \
              query_range, hmm_range, total_length \
            = line.split()
          sys.stdout.write(",%s,%s,%s," % (score, ss, cols))
          score = float(score)
          if score > max_hhalign_score:
            max_hhalign_score = score
          ss = float(ss)
          if score > max_hhalign_ss:
            max_hhalign_ss = ss
          cols = int(cols)
          if cols > max_hhalign_cols:
            max_hhalign_cols = cols
          break
    else:
      sys.stdout.write(",,,,")
    num_identities = 0
    sum_blosum62_score = 0.0
    sum_hmmsum_score = 0.0
    num_aligned_pairs = float(num_aligned_pairs)
    num_secondary_structure_agreeing = 0
    weighted_num_secondary_structure_agreeing = 0
    total_secondary_structure_weights = 0
    sum_secondary_structure_dot_product = 0
    sum_solvent_accessibility_difference = 0
    total_solvent_accessibility_weights = 0
    sum_weighted_solvent_accessibility_difference = 0
    sum_phi_difference = 0.0
    sum_psi_difference = 0.0
    for i, j in aln.pairsNoGaps:
      if seqX[i] == seqY[j]:
        num_identities += 1
      try:
        sum_blosum62_score += blosum62[(seqX[i], seqY[j])]
      except KeyError:
        sum_blosum62_score += blosum62[(seqY[j], seqX[i])]
      sum_hmmsum_score += hmmsum_score(seqX[i], seqY[j], 
                                      x_gammas[i], y_gammas[j], hmmsum_matrices)
      ss_weight = ss_X.conf[i] * ss_Y.conf[j]
      total_secondary_structure_weights += ss_weight
      if ss_X.ss[i] == ss_Y.ss[j]:
        num_secondary_structure_agreeing += 1
        weighted_num_secondary_structure_agreeing += ss_weight
      sum_secondary_structure_dot_product += \
        sum(ss_X.prob[i][k] * ss_Y.prob[j][k] for k in range(3))
      solvent_accessibility_difference = abs(sa_X[i] - sa_Y[j])
      sum_solvent_accessibility_difference += solvent_accessibility_difference
      sa_weight = sa_conf_X[i] * sa_conf_Y[j]
      total_solvent_accessibility_weights += sa_weight
      sum_weighted_solvent_accessibility_difference \
          += sa_weight * solvent_accessibility_difference
      phi_difference = abs(phi_X[i] - phi_Y[j])
      phi_difference = min(phi_difference, 360.0 - phi_difference)
      sum_phi_difference += phi_difference
      psi_difference = abs(psi_X[i] - psi_Y[j])
      psi_difference = min(psi_difference, 360.0 - psi_difference)
      sum_psi_difference += psi_difference
    if num_aligned_pairs == 0:
      percent_id = 0.0
      mean_blosum62_score = 0.0
      mean_hmmsum_score = 0.0
      fraction_ss_agreement = 0.0
      mean_secondary_structure_dot_product = 0.0
      mean_solvent_accessibility_difference = 100.0
      mean_weighted_solvent_accessibility_difference = 100.0
      mean_phi_difference = 180.0
      mean_psi_difference = 180.0
    else:
      percent_id = num_identities / num_aligned_pairs
      mean_blosum62_score = sum_blosum62_score / num_aligned_pairs
      mean_hmmsum_score = sum_hmmsum_score / num_aligned_pairs
      fraction_ss_agreement = num_secondary_structure_agreeing / \
                              num_aligned_pairs
      mean_secondary_structure_dot_product = \
          sum_secondary_structure_dot_product / num_aligned_pairs
      mean_solvent_accessibility_difference = \
          sum_solvent_accessibility_difference / \
            num_aligned_pairs
      mean_phi_difference = sum_phi_difference / num_aligned_pairs
      mean_psi_difference = sum_psi_difference / num_aligned_pairs
    if total_secondary_structure_weights == 0:
      weighted_fraction_ss_agreement = 0.0
    else:
      weighted_fraction_ss_agreement = \
          float(weighted_num_secondary_structure_agreeing) / \
          total_secondary_structure_weights
    if total_solvent_accessibility_weights == 0:
      mean_weighted_solvent_accessibility_difference = 10000
    else:
      mean_weighted_solvent_accessibility_difference = \
          float(sum_weighted_solvent_accessibility_difference) / \
              total_solvent_accessibility_weights
    if percent_id > max_percent_id:
      max_percent_id = percent_id
    if mean_blosum62_score > max_mean_blosum62_score:
      max_mean_blosum62_score = mean_blosum62_score
    if mean_hmmsum_score > max_mean_hmmsum_score:
      max_mean_hmmsum_score = mean_hmmsum_score
    if fraction_ss_agreement > max_fraction_ss_agreement:
      max_fraction_ss_agreement = fraction_ss_agreement
    if weighted_fraction_ss_agreement > max_weighted_fraction_ss_agreement:
      max_weighted_fraction_ss_agreement = weighted_fraction_ss_agreement
    if total_secondary_structure_weights > \
        max_total_secondary_structure_weights:
      max_total_secondary_structure_weights = total_secondary_structure_weights
    if mean_secondary_structure_dot_product > max_mean_ss_dot_product:
      max_mean_ss_dot_product = mean_secondary_structure_dot_product
    if mean_solvent_accessibility_difference < \
        min_mean_solvent_accessibility_difference:
      min_mean_solvent_accessibility_difference \
          = mean_solvent_accessibility_difference
    if mean_weighted_solvent_accessibility_difference < \
        min_weighted_mean_solvent_accessibility_difference:
      min_weighted_mean_solvent_accessibility_difference \
          = mean_weighted_solvent_accessibility_difference
    if total_solvent_accessibility_weights > \
        max_total_solvent_accessibility_weights:
      max_total_solvent_accessibility_weights \
          = total_solvent_accessibility_weights
    if mean_phi_difference < min_mean_phi_difference:
      min_mean_phi_difference = mean_phi_difference
    if mean_psi_difference < min_mean_psi_difference:
      min_mean_psi_difference = mean_psi_difference
    sys.stdout.write("%0.2f," % percent_id)
    sys.stdout.write("%0.3f," % mean_blosum62_score)
    sys.stdout.write("%0.3f," % mean_hmmsum_score)
    sys.stdout.write("%g," % fraction_ss_agreement)
    sys.stdout.write("%0.1f," % weighted_fraction_ss_agreement)
    sys.stdout.write("%d," % total_secondary_structure_weights)
    sys.stdout.write("%g," % mean_secondary_structure_dot_product)
    sys.stdout.write("%g," % mean_solvent_accessibility_difference)
    sys.stdout.write("%g," % mean_weighted_solvent_accessibility_difference)
    sys.stdout.write("%d," % total_solvent_accessibility_weights)
    sys.stdout.write("%0.1f," % mean_phi_difference)
    # Last feature - no comma
    sys.stdout.write("%0.1f" % mean_psi_difference)
    sys.stdout.write("\n")
  sys.stdout.write("Maximum,%g," % max_q_combined)
  sys.stdout.write("%g,%g,%g," % (max_q_developer, max_q_modeler, 
                                  max_cline_shift))
  sys.stdout.write("%d," % max_num_aligned_pairs)
  sys.stdout.write("%g," % max_yl_score)
  sys.stdout.write("%g,%g,%d," % (max_hhalign_score, max_hhalign_ss,
                                  max_hhalign_cols))
  sys.stdout.write("%g," % max_percent_id)
  sys.stdout.write("%g," % max_mean_blosum62_score)
  sys.stdout.write("%g," % max_mean_hmmsum_score)
  sys.stdout.write("%g," % max_fraction_ss_agreement)
  sys.stdout.write("%0.1f," % max_weighted_fraction_ss_agreement)
  sys.stdout.write("%d," % max_total_secondary_structure_weights)
  sys.stdout.write("%g," % max_mean_ss_dot_product)
  sys.stdout.write("%g," % min_mean_solvent_accessibility_difference)
  sys.stdout.write("%g," % min_weighted_mean_solvent_accessibility_difference)
  sys.stdout.write("%g," % max_total_solvent_accessibility_weights)
  sys.stdout.write("%g," % min_mean_phi_difference)
  # Last feature - no comma
  sys.stdout.write("%g" % min_mean_psi_difference)
  sys.stdout.write("\n")
  os.unlink(temp_afa_name)
  os.rmdir(qscore_directory)

if __name__ == '__main__':
  main()
