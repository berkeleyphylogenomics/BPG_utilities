#!/usr/bin/env python

import sys
from matchmaker.shmm_shmm_lib import *

def main():
  if len(sys.argv) < 1:
    print "Usage: %s pair" % sys.argv[0]
    sys.exit(1)

  pair = sys.argv[1]
  f = open(os.path.join(results_dir(), "%s_align" % pair,
            "%s_all_signals.csv" % pair))
  lines = [line.rstrip() for line in f.readlines()]
  f.close()
  max_line = lines[-1]
  max, max_q_combined, max_q_developer, max_q_modeler, max_cline_shift, \
  max_num_aligned_pairs, max_yl_score, max_hhalign_score, max_hhalign_ss, \
  max_hhalign_cols, max_percent_id, max_mean_blosum62_score, \
  max_mean_hmmsum_score, max_fraction_ss_agreement, \
  max_weighted_fraction_ss_agreement, max_total_secondary_structure_weights, \
  max_mean_ss_dot_product, min_mean_solvent_accessibility_difference, \
  min_weighted_mean_solvent_accessibility_difference, \
  max_total_solvent_accessibility_weights, min_mean_phi_difference, \
  min_mean_psi_difference = max_line.split(',')
  max_q_combined = float(max_q_combined)
  max_q_developer = float(max_q_developer)
  max_q_modeler = float(max_q_modeler)
  max_cline_shift = float(max_cline_shift)
  max_num_aligned_pairs = int(max_num_aligned_pairs)
  max_yl_score = float(max_yl_score)
  max_hhalign_score = float(max_hhalign_score)
  max_hhalign_ss = float(max_hhalign_ss)
  max_hhalign_cols = int(max_hhalign_cols)
  max_percent_id = float(max_percent_id)
  max_mean_blosum62_score = float(max_mean_blosum62_score)
  max_mean_hmmsum_score = float(max_mean_hmmsum_score)
  max_fraction_ss_agreement = float(max_fraction_ss_agreement)
  max_weighted_fraction_ss_agreement \
      = float(max_weighted_fraction_ss_agreement)
  max_total_secondary_structure_weights \
      = int(max_total_secondary_structure_weights)
  max_mean_ss_dot_product = float(max_mean_ss_dot_product)
  min_mean_solvent_accessibility_difference \
      = float(min_mean_solvent_accessibility_difference)
  min_weighted_mean_solvent_accessibility_difference \
      = float(min_weighted_mean_solvent_accessibility_difference)
  max_total_solvent_accessibility_weights \
      = int(max_total_solvent_accessibility_weights)
  min_mean_phi_difference = float(min_mean_phi_difference)
  min_mean_psi_difference = float(min_mean_psi_difference)
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
  for line in lines[1:-1]:
    alignment, q_combined, q_developer, q_modeler, cline_shift, \
    num_aligned_pairs, yl_score, hhalign_score, hhalign_ss, \
    hhalign_cols, percent_id, mean_blosum62_score, \
    mean_hmmsum_score, fraction_ss_agreement, \
    weighted_fraction_ss_agreement, total_secondary_structure_weights, \
    mean_ss_dot_product, mean_solvent_accessibility_difference, \
    weighted_mean_solvent_accessibility_difference, \
    total_solvent_accessibility_weights, mean_phi_difference, \
    mean_psi_difference = line.split(',')
    q_combined = float(q_combined) - max_q_combined
    q_developer = float(q_developer) - max_q_developer
    q_modeler = float(q_modeler) - max_q_modeler
    cline_shift = float(cline_shift) - max_cline_shift
    num_aligned_pairs = int(num_aligned_pairs) - max_num_aligned_pairs
    try:
      yl_score = "%g" % (float(yl_score) - max_yl_score)
    except ValueError:
      yl_score = ""
    try:
      hhalign_score = "%g" % (float(hhalign_score) - max_hhalign_score)
      hhalign_ss = "%g" % (float(hhalign_ss) - max_hhalign_ss)
      hhalign_cols = "%d" % (int(hhalign_cols) - max_hhalign_cols)
    except ValueError:
      hhalign_score = ""
      hhalign_ss = ""
      hhalign_cols = ""
    percent_id = float(percent_id) - max_percent_id
    mean_blosum62_score = float(mean_blosum62_score) - max_mean_blosum62_score
    mean_hmmsum_score = float(mean_hmmsum_score) - max_mean_hmmsum_score
    fraction_ss_agreement = float(fraction_ss_agreement) \
                              - max_fraction_ss_agreement
    weighted_fraction_ss_agreement \
        = float(weighted_fraction_ss_agreement) \
          - max_weighted_fraction_ss_agreement
    total_secondary_structure_weights \
        = int(total_secondary_structure_weights) \
          - max_total_secondary_structure_weights
    mean_ss_dot_product = float(mean_ss_dot_product) - max_mean_ss_dot_product
    mean_solvent_accessibility_difference \
        = min_mean_solvent_accessibility_difference \
          - float(mean_solvent_accessibility_difference) 
    weighted_mean_solvent_accessibility_difference \
        = min_weighted_mean_solvent_accessibility_difference \
          - float(weighted_mean_solvent_accessibility_difference)
    total_solvent_accessibility_weights \
        = int(total_solvent_accessibility_weights) \
            - max_total_solvent_accessibility_weights
    mean_phi_difference = min_mean_phi_difference - float(mean_phi_difference)
    mean_psi_difference = min_mean_psi_difference - float(mean_psi_difference)
    sys.stdout.write("%s," % alignment)
    sys.stdout.write("%g," % q_combined)
    sys.stdout.write("%g,%g,%g," % (q_developer, q_modeler, cline_shift))
    sys.stdout.write("%d," % num_aligned_pairs)
    sys.stdout.write("%s," % yl_score) 
    sys.stdout.write("%s,%s,%s," % (hhalign_score, hhalign_ss, hhalign_cols))
    sys.stdout.write("%0.2f," % percent_id)
    sys.stdout.write("%0.3f," % mean_blosum62_score)
    sys.stdout.write("%0.3f," % mean_hmmsum_score)
    sys.stdout.write("%g," % fraction_ss_agreement)
    sys.stdout.write("%0.1f," % weighted_fraction_ss_agreement)
    sys.stdout.write("%d," % total_secondary_structure_weights)
    sys.stdout.write("%g," % mean_ss_dot_product)
    sys.stdout.write("%g," % mean_solvent_accessibility_difference)
    sys.stdout.write("%g," % weighted_mean_solvent_accessibility_difference)
    sys.stdout.write("%d," % total_solvent_accessibility_weights)
    sys.stdout.write("%0.1f," % mean_phi_difference)
    # Last feature - no comma
    sys.stdout.write("%0.1f" % mean_psi_difference)
    sys.stdout.write("\n")

if __name__ == '__main__':
  main()
