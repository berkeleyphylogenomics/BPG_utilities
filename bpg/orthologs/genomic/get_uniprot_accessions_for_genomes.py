#!/usr/bin/env python

import os, sys, glob

def main():
  indir = '/clusterfs/ohana/external/genomes/QuestForOrthologs/Release5'
  os.chdir(indir)
  taxon_ids = [
     13684, # Phaeosphaeria nodorum
      1148, # Synechocystis sp. (strain ATCC 27184 / PCC 6803 / N-1)
     10116, # Rattus norvegicus
     10090, # Mus musculus
      1299, # Deinococcus radiodurans
      1423, # Bacillus subtilis
    145481, # Physcomitrella patens subsp. patens
       173, # Leptospira interrogans
      1773, # Mycobacterium tuberculosis
    184922, # Giardia intestinalis (strain ATCC 50803 / WB clone C6)
      1902, # Streptomyces coelicolor
      2190, # Methanocaldococcus jannaschii
      2214, # Methanosarcina acetivorans
      2287, # Sulfolobus solfataricus
      2336, # Thermotoga maritima
    265606, # Rhodopirellula baltica
    289376, # Thermodesulfovibrio yellowstonii (strain ATCC 51303 / DSM 11347 / YP87)
    315277, # Chlamydia trachomatis serovar A (strain HAR-13 / ATCC VR-571B)
    324602, # Chloroflexus aurantiacus (strain ATCC 29366 / DSM 635 / J-10-fl)
     33072, # Gloeobacter violaceus
     35128, # Thalassiosira pseudonana
     35554, # Geobacter sulfurreducens
     36329, # Plasmodium falciparum (isolate 3D7)
      3702, # Arabidopsis thaliana
    374847, # Korarchaeum cryptofilum (strain OPF8)
       375, # Bradyrhizobium japonicum
    381754, # Pseudomonas aeruginosa (strain PA7)
     44689, # Dictyostelium discoideum
    451804, # Aspergillus fumigatus (strain CEA10 / CBS 144.89 / FGSC A1163)
     45351, # Nematostella vectensis
    478009, # Halobacterium salinarum (strain ATCC 29341 / DSM 671 / R1)
      4896, # Schizosaccharomyces pombe
      4932, # Saccharomyces cerevisiae
      4952, # Yarrowia lipolytica
      5141, # Neurospora crassa
    515635, # Dictyoglomus turgidum (strain Z-1310 / DSM 6724)
      5207, # Cryptococcus neoformans
      5270, # Ustilago maydis
      5476, # Candida albicans
      5664, # Leishmania major
      5722, # Trichomonas vaginalis
      6183, # Schistosoma mansoni
      6239, # Caenorhabditis elegans
     63363, # Aquifex aeolicus
    665079, # Sclerotinia sclerotiorum (strain ATCC 18683 / 1980 / Ss-1)
     69014, # Thermococcus kodakarensis KOD1
      6945, # Ixodes scapularis
      7165, # Anopheles gambiae
      7227, # Drosophila melanogaster
     76856, # Fusobacterium nucleatum subsp. nucleatum
      7739, # Branchiostoma floridae
      7955, # Danio rerio
     81824, # Monosiga brevicollis
       818, # Bacteroides thetaiotaomicron
     83333, # Escherichia coli (strain K12)
      9031, # Gallus gallus
      9606, # Homo sapiens
  ]
  for taxon_id in taxon_ids:
    file = glob.glob('%d_*.fasta' % taxon_id)[0]
    f = open(file)
    for line in f.readlines():
      if line[0] == '>':
        identifier = line[1:].strip().split()[0]
        if identifier[0:9] != 'UniProtKB':
          sys.stderr.write("Unrecognized identifier in %d: %s\n" 
                            % (taxon_id, line.strip()))
        else:
          uniprot_accession = identifier.split(':')[1]
          sys.stdout.write("%d,%s\n" % (taxon_id, uniprot_accession))

if __name__ == '__main__':
  main()
