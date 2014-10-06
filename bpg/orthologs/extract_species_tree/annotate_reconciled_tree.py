#!/usr/bin/python

import sys, string, re, MySQLdb
from connect_db import *

db = MySQLdb.connect( db=phylofacts_db_name, host=phylofacts_db_host, \
                      user=phylofacts_db_user, passwd=phylofacts_db_passwd )
cur = db.cursor( MySQLdb.cursors.DictCursor )

def annotate_nhx_element(bpg_accession, element):
  tags = element.split(':')
  if len(tags) >= 1:
    dbrefs = tags[0].split('_bpgseq')
    if len(dbrefs) > 1:
      taxon_id = dbrefs[0]
      sequence_id = dbrefs[1]
      sql = """ SELECT  scientific_name, common_name, 
                        subfamily_node, subfamily_alignment_id,
                        uniprot_accession, uniprot_de
                FROM    
                        protein_sequence, 
                        alignment_protein_sequence, 
                        alignment,
                        book, 
                        genbank_uniprot, 
                        ncbi_taxonomy
                WHERE   protein_sequence.id = protein_sequence_id
                AND     alignment_id = alignment.id
                AND     book_id = book.id
                AND     scopid_or_bpgid LIKE '%s%%'
                AND     subfamily_alignment_id = 0
                AND     genbank_uniprot_id = genbank_uniprot.id
                AND     ncbi_taxonomy.id = %s
                AND     sequence_id = %s
                ORDER BY subfamily_alignment_id DESC
            """ % (bpg_accession, taxon_id, sequence_id)
      num_rows = cur.execute(sql)
      if num_rows > 0:
        row = cur.fetchone()
        description_elements = []
        if row['subfamily_alignment_id'] != 0 and row['subfamily_node'] != None:
          description_elements.append(row['subfamily_node'])
        if row['uniprot_accession'] != None:
          description_elements.append(row['uniprot_accession'])
        if row['uniprot_de'] != None:
          description_elements.append(string.replace(row['uniprot_de'], ' ',
                                                                        '_'))
        description = '_'.join(description_elements)
        if 'common_name' in row and row['common_name'] != None:
          taxon = "%s_(%s)" % (row['scientific_name'], row['common_name'])
        else:
          taxon = row['scientific_name']
        taxon = string.replace(taxon, ' ', '_')
        return ':'.join(['_'.join([taxon, description]), ':'.join(tags[1:])])
  return element

def main():
  if len(sys.argv) < 3:
    print ("Usage: %s bpg_accession reconciled_tree_file") % sys.argv[0]
    sys.exit(0)

  bpg_accession = sys.argv[1]

  f = open(sys.argv[2])
  nhx_tree = f.read()
  f.close()


  nhx_element_re = re.compile('([(,])')
  print ''.join([annotate_nhx_element(bpg_accession, element)
                for element in nhx_element_re.split(nhx_tree)])


if __name__ == '__main__':
  main()
