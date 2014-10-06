#!/usr/bin/env python2.6

from cogent.core.tree import PhyloNode
from cogent import LoadTree
from cogent.parse.tree import DndParser

def main():

  f = open('/clusterfs/ohana/external/SILVA/LTP_release_104/LTPs104_SSU_tree.newick')
  tree_string = f.read()
  f.close()

  tree = DndParser(tree_string, PhyloNode)
  taxon_id_of_name = {}
  taxon_id_of_name['Deinococcus_radiodurans__Y11332__Deinococcaceae'] = 1299
  taxon_id_of_name['Bacillus_subtilis_subsp._subtilis__AJ276351__Bacillaceae'] = 1423
  taxon_id_of_name['Leptospira_interrogans__Z12817__Leptospiraceae'] = 173
  taxon_id_of_name['Mycobacterium_tuberculosis__X58890__Mycobacteriaceae'] = 1773
  taxon_id_of_name['Streptomyces_coelicoflavus__AB184650__Streptomycetaceae'] = 1902
  taxon_id_of_name['Methanocaldococcus_jannaschii___L77117__Methanocaldococcaceae'] = 2190
  taxon_id_of_name['Methanosarcina_acetivorans__M59137__Methanosarcinaceae'] = 2214
  taxon_id_of_name['Sulfolobus_solfataricus__D26490__Sulfolobaceae'] = 2287
  taxon_id_of_name['Thermotoga_maritima__M21774__Thermotogaceae'] = 2336
  taxon_id_of_name['Rhodopirellula_baltica__BX294149__Planctomycetaceae'] = 265606
  taxon_id_of_name['Thermodesulfovibrio_yellowstonii__AB231858__Nitrospiraceae'] = 289376
  taxon_id_of_name['Chlamydia_trachomatis__D89067__Chlamydiaceae'] = 315277
  taxon_id_of_name['Chloroflexus_aurantiacus__D38365__Chloroflexaceae'] = 324602
  taxon_id_of_name['Geobacter_sulfurreducens__U13928__Geobacteraceae'] = 35554
  taxon_id_of_name['Bradyrhizobium_japonicum__U69638__Bradyrhizobiaceae'] = 375
  taxon_id_of_name['Pseudomonas_aeruginosa__X06684__Pseudomonadaceae'] = 381754
  taxon_id_of_name['Halobacterium_salinarum__AJ496185__Halobacteriaceae'] = 478009
  taxon_id_of_name['Dictyoglomus_turgidum__CP001251__Dictyoglomaceae'] = 515635
  taxon_id_of_name['Aquifex_pyrophilus__M83548__Aquificaceae'] = 63363
  taxon_id_of_name['Thermococcus_kodakarensis__D38650__Thermococcaceae'] = 69014
  taxon_id_of_name['Fusobacterium_nucleatum_subsp._nucleatum__AE009951__Fusobacteriaceae'] = 76856
  taxon_id_of_name['Bacteroides_thetaiotaomicron___AE015928__Bacteroidaceae'] = 818
  taxon_id_of_name['Escherichia_coli__X80725__Enterobacteriaceae'] = 83333
  node_of_taxon_id = {}
  for name in taxon_id_of_name:
    node_of_taxon_id[taxon_id_of_name[name]] = tree.getNodeMatchingName(name)
  max_distance = 0.0
  for taxon_id1 in node_of_taxon_id:
    for taxon_id2 in node_of_taxon_id:
      if taxon_id1 < taxon_id2:
        distance \
            = node_of_taxon_id[taxon_id1].distance(node_of_taxon_id[taxon_id2])
        if distance > max_distance:
          max_distance = distance
        print "dist[%d][%d] = %g" % (taxon_id1, taxon_id2, distance)
  print "Maximum distance: %g" % max_distance
  scale = round(2.5 / max_distance)
  print "Scale:", scale
  for taxon_id1 in node_of_taxon_id:
    for taxon_id2 in node_of_taxon_id:
      if taxon_id1 < taxon_id2:
        threshold \
            = node_of_taxon_id[taxon_id1].distance(node_of_taxon_id[taxon_id2])
        threshold *= scale
        threshold = round(threshold * 8)
        threshold = 0.125 * threshold
        print "threshold_of_taxon_pair[%d][%d] = %g" \
          % (taxon_id1, taxon_id2, threshold)

  print "Threholds with Eukaryotes"
  for taxon_id1 in [3702, 4896, 4932, 6239, 7227, 7955, 9606, 10090, 44689]:
    for taxon_id2 in node_of_taxon_id:
      if taxon_id1 < taxon_id2:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id1, taxon_id2)
      else:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id2, taxon_id1)

  for taxon_id1 in [3702, 4896, 4932, 6239, 7227, 7955, 9606, 10090, 44689]:
    for taxon_id2 in [1148, 33072, 374847]:
      if taxon_id1 < taxon_id2:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id1, taxon_id2)
      else:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id2, taxon_id1)

  print "Thresholds with more Eukaryotes"
  for taxon_id1 in node_of_taxon_id:
    for taxon_id2 in [10116, 9031, 81824, 7739, 7165, 6945, 665079, 6183,
    5476, 5722, 5664, 5270, 5207, 5141, 4952, 45351, 451804, 36329,
    35128, 184922, 145481, 13684]:
      if taxon_id1 < taxon_id2:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id1, taxon_id2)
      else:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id2, taxon_id1)

  for taxon_id1 in [1148, 33072, 374847]:
    for taxon_id2 in [10116, 9031, 81824, 7739, 7165, 6945, 665079, 6183,
    5476, 5722, 5664, 5270, 5207, 5141, 4952, 45351, 451804, 36329,
    35128, 184922, 145481, 13684]:
      if taxon_id1 < taxon_id2:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id1, taxon_id2)
      else:
        print "threshold_of_taxon_pair[%d][%d] = 2.5" % (taxon_id2, taxon_id1)

  print "Thresholds among Eukaryotes"
  for taxon_id1 in [5664, 5722, 35128, 36329, 184922]:
    for taxon_id2 in [3702, 4896, 4932,
    6239, 7227, 7955, 9606, 10090, 44689, 10116, 9031, 81824, 7739, 7165,
    6945, 665079, 6183, 5476, 5722, 5664, 5270, 5207, 5141, 4952, 45351,
    451804, 36329, 35128, 184922, 145481, 13684]:
      if taxon_id1 < taxon_id2:
        print "threshold_of_taxon_pair[%d][%d] = 1.5" % (taxon_id1, taxon_id2)
      elif taxon_id2 < taxon_id1:
        print "threshold_of_taxon_pair[%d][%d] = 1.5" % (taxon_id2, taxon_id1)

  taxon_id1 = 145481

  for taxon_id2 in [5664, 5722, 35128, 36329, 184922, 4896, 4932,
  6239, 7227, 7955, 9606, 10090, 44689, 10116, 9031, 81824, 7739, 7165,
  6945, 665079, 6183, 5476, 5722, 5664, 5270, 5207, 5141, 4952, 45351,
  451804, 36329, 35128, 184922, 145481, 13684]:
    if taxon_id1 < taxon_id2:
      print "threshold_of_taxon_pair[%d][%d] = 1.5" % (taxon_id1, taxon_id2)
    elif taxon_id2 < taxon_id1:
      print "threshold_of_taxon_pair[%d][%d] = 1.5" % (taxon_id2, taxon_id1)

  print "Thresholds between Fungi and non-fungi Eukaryotes"
  for taxon_id1 in [665079, 4952, 5141, 5207, 451804, 5270, 5476, 13684,
                    81824]:
    for taxon_id2 in [9606, 10090, 7955, 7227, 6239, 44689, 3702]:
      if taxon_id1 < taxon_id2:
        print "threshold_of_taxon_pair[%d][%d] = 1.5" % (taxon_id1, taxon_id2)
      else:
        print "threshold_of_taxon_pair[%d][%d] = 1.5" % (taxon_id2, taxon_id1)
if __name__ == '__main__':
  main()
