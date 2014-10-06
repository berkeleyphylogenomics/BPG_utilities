#!/usr/bin/env python

def main():

  # RSD 2010/10/05: I ran the following three commands:
  # echo "COPY (SELECT phogt_loose_node_id, uniprot_identifier FROM tree_node, sequence_header, uniprot WHERE sequence_header_id = sequence_header.id AND uniprot_id = uniprot.id AND uniprot.taxon_id = 9606 ORDER BY phogt_loose_node_id) TO STDOUT WITH CSV HEADER" | psql -o /home/ruchira/human_loose_phogs.csv -h db pfacts003_test -U bpg_user
  # echo "COPY (SELECT phogt_loose_node_id, uniprot_identifier FROM tree_node, sequence_header, uniprot WHERE sequence_header_id = sequence_header.id AND uniprot_id = uniprot.id AND uniprot.taxon_id = 7227 ORDER BY phogt_loose_node_id) TO STDOUT WITH CSV HEADER" | psql -o /home/ruchira/fruitfly_loose_phogs.csv -h db pfacts003_test -U bpg_user
  # cat human_loose_phogs.csv fruitfly_loose_phogs.csv | sort -t ',' -k 1,2 > human_fruitfly_loose_phogs.csv
  # This produced a list of all loose PHOGs containing human or Drosophila
  # sequences.  Here is a segment of human_fruitfly_loose_phogs.csv:
  # 10006563,Q4G0X3_HUMAN
  # 10006563,Q4VXY0_HUMAN
  # 10006563,Q8IQN2_DROME
  # 10006563,Q8MQK2_DROME
  
  # Now we parse this file to produce PHOGs that contain both human and
  # Drosophila sequences.

  f = open("/home/ruchira/human_fruitfly_loose_phogs.csv")
  lines = f.readlines()
  f.close()

  orthologs_in_phog = {}

  for line in lines:
    phog, ortholog = line.rstrip().split(',')
    species = ortholog.split('_')[1]
    if phog not in orthologs_in_phog:
      orthologs_in_phog[phog] = {}
    if species not in orthologs_in_phog[phog]:
      orthologs_in_phog[phog][species] = set()
    orthologs_in_phog[phog][species].add(ortholog)

  f = open('/clusterfs/ohana/sandbox/human_fruitfly_loose_orthologs.csv','w')
  f.write('"H. sapien co-orthologs","D. melanogaster co-orthologs"\n')
  for phog in orthologs_in_phog.keys():
    if 'HUMAN' in orthologs_in_phog[phog] and \
        'DROME' in orthologs_in_phog[phog]:
      f.write('"%s","%s"\n' % (';'.join(list(orthologs_in_phog[phog]['HUMAN'])),
                          ';'.join(list(orthologs_in_phog[phog]['DROME']))))
  f.close()

if __name__ == '__main__':
  main()
