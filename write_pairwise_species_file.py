#!/usr/bin/env python

from optparse import OptionParser
import os, glob, re, sys, math, subprocess, cPickle

def write_protein_preamble(taxon_id, taxon_name, outf):
  outf.write('  <species name="%s" NCBITaxId="%d">\n' % (taxon_name, taxon_id))
  outf.write('    <database name="UniProt" version="2011_04" \n')
  outf.write('      protLink="http://www.uniprot.org/uniprot/">\n')
  outf.write('      <genes>\n')

def write_protein_postamble(outf):
  outf.write('      </genes>\n')
  outf.write('    </database>\n')
  outf.write('  </species>\n')

def write_ortholog_preamble(outf):
  outf.write('  <scores>\n')
  outf.write('    <scoreDef id="PHOG" ')
  outf.write(
    'desc="exp(-(minimum threshold of smallest containing PHOG))" />\n')
  outf.write('  </scores>\n')
  outf.write('  <groups>\n')

def write_ortholog_postamble(outf):
  outf.write('  </groups>\n')
  outf.write('</orthoXML>\n')

def write_orthologs(taxon_id1, taxon_id2, threshold, species_1_f, species_2_f,
                    ortholog_f):
  family_of_tree = {}
  f = open('/home/ruchira/good_family_tree_ids_20110608.txt')
  f.readline()
  for line in f.readlines():
    tree_id, family_accession = line.strip().split(',')
    family_of_tree[tree_id] = family_accession
  f.close()
  taxon_id1_str = str(taxon_id1)
  taxon_id2_str = str(taxon_id2)
  query_taxon_re = re.compile('[\(,][A-Z0-9]+:(\d+)\)')
  identifier_re = re.compile('(\d+) \(([A-Z0-9]+):(\d+)[,A-Z0-9:]*\)')
  ortholog_sep_re = re.compile('};')
  support_sep_re = re.compile(' <= ')
  support_re = re.compile('\(PHOG(\d+)_(\d+), ([0-9]+\.?[0-9e\-]*)\),')
  written_query_genes = set()
  written_target_genes = set()
  num_processed_lines = 0
  def process_line(line, target_taxon_str, query_species_f, target_species_f):
    line = line.strip()
    m = identifier_re.match(line)
    query_uniprot_id = m.group(1)
    query_uniprot_accession = m.group(2)
    query_taxon = m.group(3)
    # Skip past match and ': '
    line = line[len(m.group(0))+2:]
    ortholog_specs = ortholog_sep_re.split(line)
    foundOrthologs = False
    supporting_phogs = {}
    accession_of_id = {}
    for ortholog_spec in ortholog_specs:
      if len(ortholog_spec[1:]) == 0:
        continue
      # Skip past the initial '{'
      identifier, support = support_sep_re.split(ortholog_spec[1:])
      om = identifier_re.match(identifier)
      if om and om.group(3) == target_taxon_str:
        ortholog_id = om.group(1)
        accession_of_id[ortholog_id] = om.group(2)
        supporting_phogs[ortholog_id] = set()
        support_match = support_re.match(support)
        while support_match and len(support_match.group(0)) > 0:
          min_threshold = float(support_match.group(3))
          if min_threshold <= threshold:
            tree_id = support_match.group(1).lstrip('0')
            left_id = support_match.group(2)
            family_id = family_of_tree[tree_id]
            supporting_phogs[ortholog_id].add((min_threshold, family_id, 
                                              tree_id, left_id))
          support = support[len(support_match.group(0)):]
          support_match = support_re.match(support)
        if len(supporting_phogs[ortholog_id]) > 0:
          supporting_phogs[ortholog_id] = list(supporting_phogs[ortholog_id])
          supporting_phogs[ortholog_id].sort()
          foundOrthologs = True
      # Otherwise this ortholog is not of interest to us
    if foundOrthologs:
      for ortholog_id in supporting_phogs:
        if len(supporting_phogs[ortholog_id]) > 0:
          if query_uniprot_id not in written_query_genes:
            query_species_f.write('        <gene id="%s" protId="%s" />\n' 
                                  % (query_uniprot_id, query_uniprot_accession))
            written_query_genes.add(query_uniprot_id)
          ortholog_f.write('    <orthologGroup id="%s_%s">\n' % (query_uniprot_id, ortholog_id))
          score = math.exp(-supporting_phogs[ortholog_id][0][0])
          ortholog_f.write('      <geneRef id="%s">\n' % query_uniprot_id)
          ortholog_f.write('      </geneRef>\n')
          target_gene_id = int(ortholog_id)
          if target_gene_id not in written_target_genes:
            target_species_f.write('        <gene id="%s" protId="%s" />\n'
                              % (ortholog_id, accession_of_id[ortholog_id]))
            written_target_genes.add(target_gene_id)
          ortholog_f.write('      <geneRef id="%s">\n' % ortholog_id)
          ortholog_f.write('        <score id="PHOG" value="%0.4f" />\n' 
                            % score)
          ortholog_f.write('        <notes>\n')
          for (min_threshold, family_id, tree_id, left_id) \
              in supporting_phogs[ortholog_id]:
            ortholog_f.write('          <supporting_phog ')
            ortholog_f.write('family="%s" ' % family_id)
            ortholog_f.write('phog="PHOG%s_%s" ' % (tree_id, left_id))
            ortholog_f.write('minimum_threshold="%0.4f" ' % min_threshold)
            ortholog_f.write('/>\n')
          ortholog_f.write('        </notes>\n')
          ortholog_f.write('      </geneRef>\n')
          ortholog_f.write('    </orthologGroup>\n')

  ortholog_files = glob.glob('/clusterfs/vasudha/bpg/OrthologsForQuest/*_of_37')
  print "Found %d ortholog files" % len(ortholog_files)
  for file in ortholog_files:
    f = open(file)
    for line in f.readlines():
      m = query_taxon_re.search(line)
      if m is None:
        print "Unmatching line:", line
      query_protein_taxon = m.group(1)
      if query_protein_taxon == taxon_id1_str:
        process_line(line, taxon_id2_str, species_1_f, species_2_f)
        num_processed_lines += 1
    f.close()

def main():
  taxa = {
     10090: 'Mus_musculus',
      1299: 'Deinococcus_radiodurans',
      1423: 'Bacillus_subtilis',
       173: 'Leptospira_interrogans',
      1773: 'Mycobacterium_tuberculosis',
      1902: 'Streptomyces_coelicolor',
      2190: 'Methanocaldococcus_jannaschii',
      2214: 'Methanosarcina_acetivorans',
      2287: 'Sulfolobus_solfataricus',
      2336: 'Thermotoga_maritima',
    265606: 'Rhodopirellula_baltica',
    289376: 'Thermodesulfovibrio_yellowstonii',
    315277: 'Chlamydia_trachomatis',
    324602: 'Chloroflexus_aurantiacus',
     35554: 'Geobacter_sulfurreducens',
      3702: 'Arabidopsis_thaliana',
       375: 'Bradyrhizobium_japonicum',
    381754: 'Pseudomonas_aeruginosa',
     44689: 'Dictyostelium_discoideum',
    478009: 'Halobacterium_salinarum',
      4896: 'Schizosaccharomyces_pombe',
      4932: 'Saccharomyces_cerevisiae',
    515635: 'Dictyoglomus_turgidum',
      6239: 'Caenorhabditis_elegans',
     63363: 'Aquifex_aeolicus',
     69014: 'Thermococcus_kodakarensis',
      7227: 'Drosophila_melanogaster',
     76856: 'Fusobacterium_nucleatum',
      7955: 'Danio_rerio',
       818: 'Bacteroides_thetaiotaomicron',
     83333: 'Escherichia_coli',
      9606: 'Homo_sapiens',
      1148: 'Synechocystis_sp.',
    374847: 'Korarchaeum_cryptofilum',
     33072: 'Gloeobacter_violaceus',
      5664: 'Leishmania_major',
      5722: 'Trichomonas_vaginalis',
     35128: 'Thalassiosira_pseudonana',
     36329: 'Plasmodium_falciparum',
    184922: 'Giardia_intestinalis',
    145481: 'Physcomitrella_patens',
     81824: 'Monosiga_brevicollis',
     13684: 'Phaeosphaeria_nodorum',
    665079: 'Sclerotinia_sclerotiorum',
      4952: 'Yarrowia_lipolytica',
      5141: 'Neurospora_crassa',
      5207: 'Cryptococcus_neoformans',
    451804: 'Aspergillus_fumigatus',
      5270: 'Ustilago_maydis',
      5476: 'Candida_albicans',
     45351: 'Nematostella_vectensis',
      6183: 'Schistosoma_mansoni',
      6945: 'Ixodes_scapularis',
      7165: 'Anopheles_gambiae',
      7739: 'Branchiostoma_floridae',
      9031: 'Gallus_gallus',
     10116: 'Rattus_norvegicus',
        }
  usage = "%prog [options] taxon_id1 taxon_id2 threshold"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-t", "--threshold", dest="threshold", 
                        default=-1.0,
                        help="PHOG threshold to use")
  opt_parser.add_option("-x", "--xmloutput", dest="produce_xml_output",
                      action="store_true", default=True,
                      help="Whether to write the output in orthoXML format.")
  (options, args) = opt_parser.parse_args()
  if len(args) < 2:
    opt_parser.error('Incorrect number of arguments')
  def taxon_error(opt_parser):
    opt_parser.error('Taxon_id must be one of %s'
      % ','.join(['%d [%s]' % (taxon_id, taxon_name) for (taxon_id, taxon_name)
                  in taxa.items()]))
  try:
    taxon_id1 = int(args[0])
  except ValueError:
    opt_parser.error('taxon_id1 should be an NCBI taxonomy id')
  if taxon_id1 not in taxa:
    taxon_error(opt_parser)
  taxon_name1 = taxa[taxon_id1]
  try:
    taxon_id2 = int(args[1])
  except ValueError:
    opt_parser.error('taxon_id2 should be an NCBI taxonomy id')
  if taxon_id2 not in taxa:
    taxon_error(opt_parser)
  taxon_name2 = taxa[taxon_id2]
  if options.threshold < 0.0:
    f = open('/clusterfs/vasudha/bpg/OrthologsForQuest/threshold_of_taxon_pair.pkl')
    threshold_of_taxon_pair = cPickle.load(f)
    f.close()
    threshold = threshold_of_taxon_pair[min(taxon_id1,
                                taxon_id2)][max(taxon_id1, taxon_id2)]
  else:
    threshold = options.threshold
  species_1_f = open('%s_genes.xml' % taxon_name1, 'w')
  species_2_f = open('%s_genes.xml' % taxon_name2, 'w')
  ortholog_f = open('%s-%s_orthologs_piece.xml' % (taxon_name1, taxon_name2), 
                      'w')
  write_protein_preamble(taxon_id1, taxon_name1, species_1_f)
  write_protein_preamble(taxon_id2, taxon_name2, species_2_f)
  write_ortholog_preamble(ortholog_f)
  write_orthologs(taxon_id1, taxon_id2, threshold, species_1_f, species_2_f, 
                  ortholog_f)
  write_protein_postamble(species_1_f)
  species_1_f.close()
  write_protein_postamble(species_2_f)
  species_2_f.close()
  write_ortholog_postamble(ortholog_f)
  ortholog_f.close()
  outf = open("%s-%s_%g_orthologs.xml" % (taxon_name1, taxon_name2, threshold),
              "w")
  outf.write('<?xml version="1.0" encoding="utf-8" ?>\n')
  outf.write('<orthoXML xmlns="http://orthoXML.org/2011/" version = "0.3" ')
  outf.write('origin="PHOG" originVersion="2.0" ')
  outf.write('xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ')
  outf.write('xsi:schemaLocation="http://orthoXML.org/2011/ ')
  outf.write('http://www.orthoxml.org/0.3/orthoxml.xsd">\n')
  outf.write('  <notes>\n')
  outf.write('    Orthologs between NCBI taxon ids %d and %d at threshold %g\n'
            % (taxon_id1, taxon_id2, threshold))
  outf.write('    Using PhyloFacts 3.0 as of 2011/06/08\n')
  outf.write('  </notes>\n')
  outf.close()
  outf = open("%s-%s_%g_orthologs.xml" % (taxon_name1, taxon_name2, threshold),
              "a")
  p = subprocess.Popen(['cat',
            '%s_genes.xml' % taxon_name1,
            '%s_genes.xml' % taxon_name2,
            '%s-%s_orthologs_piece.xml' 
            % (taxon_name1, taxon_name2)],
            stdout = outf)
  status = os.waitpid(p.pid, 0)[1]
  outf.close()

if __name__ == '__main__':
  main()
