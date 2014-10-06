#!/usr/bin/python
'''
A listing of the various sql commands used to query the postgres database.

Please add your own snippets to this file and update the list in this comment (to avoid duplication)
PHOG_COLLECTION_SQL
PHOG_DATA_SQL
FAMILY_TAXA_SQL
SEQUENCE_COUNT_SQL
PFAM_ACCESSION_FROM_NAME
FAMILY_PFAM_DOMAIN_NAME
'''
########## Connection to server
import psycopg2
import psycopg2.extras
from pfacts003.utils.credentials import get_credentials

#Database connection globals
DB_NAME = 'pfacts003_test'
USER = 'webuser'
PWD = get_credentials(USER)
def connect_to_server():
    """
    Connects to postgres database and returns the cursor.
    """
    conn = psycopg2.connect("dbname='%s' user='%s' host='db' password='%s'" % (DB_NAME, USER, PWD))
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    return cur

def sql_results(cur, sql_query, parameter_tuple):
    """
    Run the sql query and return results as a list of dictionaries.
    """
    cur.execute(sql_query, parameter_tuple)
    return cur.fetchall()

def sql_results_without_params(cur, sql_query):
    """
    Run the sql query and return results as a list of dictionaries.
    """
    cur.execute(sql_query)
    return cur.fetchall()

##########
########## SQL patterns
#Get the phogs belonging to a family
PHOG_COLLECTION_SQL = """
SELECT family.id as bpg_accession, tree.id as tree_id, tree_node.id as tree_node_id,
tree_node.left_id as tree_node_left_id, tree_node.is_superorthologous
FROM family, tree, tree_node
WHERE family.id = %s AND -- First select only cases where family id is the chosen one
family.canonical_tree_id = tree.id AND
tree_node.tree_id = tree.id AND -- Get the tree nodes for these trees
tree_node.is_superorthologous = 't'; -- Filter for superorthologous nodes
"""

#Get the uniprot ids, sequence headers, description and taxa found within a given PHOG.
PHOG_DATA_SQL = """
SELECT
'PHOG'||lpad(CAST(tn2.tree_id as TEXT),7,'0')||'_'||lpad(CAST(tn1.left_id as TEXT),5,'0')
as "PHOG identifier", -- Get phog accessions
uniprot.uniprot_identifier as "Uniprot identifier",
uniprot.accession as "Uniprot accession",
uniprot.de as "Uniprot description",
uniprot.taxon_id as "Taxon id",
uniprot_taxonomy.scientific_name as "Taxon scientific name", -- The queries are self explanatory
sequence_header.header as "Sequence header" -- Getting the sequence header for alignment bounds
FROM tree_node as tn1, tree_node as tn2,
sequence_header, uniprot, uniprot_taxonomy -- Two instances of tree_node
WHERE tn1.id = %s AND -- tree_node id to the first instance
tn2.tree_id = %s AND -- tree id of the first instance (Filter down to required values with this)
tn2.left_id >= tn1.left_id AND -- Now use the second instance to compare against values from the first
tn2.right_id <= tn1.right_id AND
tn2.sequence_header_id = sequence_header.id AND -- Finally get sequence headers, uniprot fields
sequence_header.uniprot_id = uniprot.id AND
uniprot_taxonomy.id = uniprot.taxon_id;
"""

#Get PHOGs for singleton sequence. In case phogs are not found, get a tree id and left id for a sequence
#Input is the family and the sequence accession
SINGLETON_PHOGS = """
SELECT
lpad(CAST(tree.id as TEXT),7,'0')||'_'||lpad(CAST(tree_node.left_id as TEXT),5,'0')
FROM family, tree, tree_node, sequence_header, uniprot
WHERE family.id = %s AND
	family.canonical_tree_id = tree.id AND
	tree_node.tree_id = tree.id AND
	tree_node.sequence_header_id is not NULL AND
	tree_node.sequence_header_id = sequence_header.id AND
	sequence_header.uniprot_id = uniprot.id AND
	uniprot.accession = %s;"""

#Given a taxon id, get the families 
FAMILY_TAXA_SQL = """
SELECT 	
	DISTINCT(family.id) as family_id
FROM
	family, tree, tree_node, sequence_header
WHERE
	family.active = 't' AND
	family.canonical_tree_id = tree.id AND
	tree_node.tree_id = tree.id AND 
	tree_node.sequence_header_id is not Null AND
	tree_node.sequence_header_id = sequence_header.id AND
	sequence_header.taxon_id = %s;    
"""

#Given a family accession, get the number of sequences within
SEQUENCE_COUNT_SQL = """
    select count(*) from family, tree, tree_node, sequence_header, uniprot 
    where family.id = %s and 
    family.canonical_tree_id = tree.id and 
    tree_node.tree_id = tree.id and
    tree_node.sequence_header_id is not NULL and
    tree_node.sequence_header_id = sequence_header.id and
    sequence_header.uniprot_id = uniprot.id;
    """

#Given a family accession list, get the uniprot accessions in this family
SEQUENCE_FAMILY_SQL = """
    select distinct(uniprot.accession) from family, tree, tree_node, sequence_header, uniprot 
    where family.id in %s and 
    family.canonical_tree_id = tree.id and 
    tree_node.tree_id = tree.id and
    tree_node.sequence_header_id is not NULL and
    tree_node.sequence_header_id = sequence_header.id and
    sequence_header.uniprot_id = uniprot.id;
    """
#Get pfam acccession when the overall_pfam_version and name are provided
PFAM_ACCESSION_FROM_NAME = """
   SELECT accession FROM PFAM WHERE name=%s AND overall_pfam_version=%s"""

#For a pfam domain accession get all the pfam families using this domain as its seed.
#Renamed from FAMILY_PFAM_DOMAIN_NAME
PFAM_ACCESSION_PFAM_FAMILY = """
SELECT distinct(family.id)
from pfam, hmm, hmm_consensus, tree_node_consensus, tree_node, tree, family, family_type
where pfam.accession = %s AND
pfam.id = hmm.pfam_id AND
hmm.id = hmm_consensus.hmm_id AND
hmm_consensus.id = tree_node_consensus.id AND
tree_node_consensus.method = 'hmm' AND
tree_node_consensus.tree_node_id = tree_node.id AND
tree_node.tree_id = tree.id AND
family.canonical_tree_id = tree.id AND
family.active = 't' AND
family.family_type_id = family_type.id AND
family_type.id = 'C';"""

# for a pfam domain accession get all the GHG families using this domain as its seed.
PFAM_ACCESSION_GHG_FAMILY = """
SELECT 
	family.id, sequence_hmm.sequence_start, sequence_hmm.sequence_end
FROM 
	pfam, hmm, sequence_hmm, tree_node_consensus, tree_node, tree, family
WHERE
	pfam.accession = %s AND
	hmm.pfam_id = pfam.id AND
	sequence_hmm.hmm_id = hmm.id AND
	sequence_hmm.sequence_id = tree_node_consensus.sequence_id AND
	tree_node.id = tree_node_consensus.tree_node_id AND
	tree_node.left_id = 1 AND 
	tree.id = tree_node.tree_id AND
	tree.id = family.canonical_tree_id AND
	family.active = 't' AND
        family.family_type_id = 'G';
"""

# For a GHG family, get the pfam domains with sequence start and end values
GHG_FAMILY_PFAM_DOMAIN_NAME = """
SELECT 
	pfam.name, sequence_hmm.sequence_start, sequence_hmm.sequence_end
FROM
	family, tree, tree_node, tree_node_consensus, sequence, sequence_hmm, hmm, pfam
WHERE
	family.id = %s AND
	family.family_type_id = 'G' AND
	tree.id = family.canonical_tree_id AND
	tree_node.tree_id = tree.id AND
	tree_node.left_id = 1 AND -- Till here, get the canonical root node
	tree_node_consensus.tree_node_id = tree_node.id AND -- Get the root consensus
	tree_node_consensus.sequence_id = sequence.id AND
	sequence_hmm.sequence_id = sequence.id AND
	hmm.id = sequence_hmm.hmm_id AND
	hmm.pfam_id is not null AND
	pfam.id = hmm.pfam_id;
"""

#Get seguid for accessions from the uniprot table
ACCESSION_SEGUID = """
SELECT seguid from uniprot where accession = %s;
"""

#Get phogs from uniprot accession
PHOG_ACCESSION_SQL = """
SELECT
	tn2.id as tree_node_id, tn2.tree_id as tree_id, tn2.left_id as tree_node_left_id
FROM
	uniprot, sequence_header, tree_node as tn1,
	tree_node as tn2, family,tree 
WHERE 
	uniprot.accession = %s AND 
	sequence_header.uniprot_id = uniprot.id AND 
	tn1.sequence_header_id = sequence_header.id AND
	tn2.tree_id = tn1.tree_id AND
	tn2.left_id <= tn1.left_id AND
	tn2.right_id >= tn1.right_id AND
	tn2.is_superorthologous='t' AND
	tree.id = tn2.tree_id AND
	family.canonical_tree_id = tree.id AND
	family.active = 't';"""

# Given a sequence, get a list of families that it is in.
FAMILY_SEQUENCE_SQL = """
SELECT
    distinct(family.id)
FROM
    uniprot, family, tree, tree_node, sequence_header
WHERE
    uniprot.accession = %s AND
    sequence_header.uniprot_id = uniprot.id AND
    tree_node.sequence_header_id = sequence_header.id AND
    tree_node.sequence_header_id is not NULL AND
    tree.id = tree_node.tree_id AND
    family.canonical_tree_id = tree.id AND
    family.status != 'bad' AND
    family.active = true;
    """

#Get go annotations from a list of sequences
SEQUENCE_LIST_GO = """
SELECT
    go_term.name, go_term.term_type, go_term.acc, go_evidence.code 
FROM
uniprot, uniprot_go, go_term, go_evidence 
WHERE
go_term.id = uniprot_go.go_term_id and
go_evidence.id = uniprot_go.go_evidence_id and
uniprot_go.uniprot_id = uniprot.id AND
uniprot.accession in %s
order by
term_type"""


#Get EC information from a list of sequences
SEQUENCE_LIST_EC = """
select
ec.id, ec.class_number, ec.subclass_number, ec.subsubclass_number,
ec.enzyme_number, ec.description
from
uniprot, uniprot_ec, ec
where
uniprot.accession in %s AND
uniprot_ec.uniprot_id = uniprot.id AND
ec.id = uniprot_ec.ec_id;
"""

#Get KEGG information given EC id
EC_LIST_KEGG = """
select
kegg_map.id, kegg_map.title
from
kegg_map_ec, kegg_map
where
kegg_map_ec.ec_id in %s AND
kegg_map_ec.kegg_map_id = kegg_map.id
"""

# Return the best GHG book for a uniprot accession
# We pick the best book based on the one that has the largest number of taxa
BEST_GHG_BOOK = """
(
SELECT
    family.id, COUNT( DISTINCT all_uniprot.taxon_id) as number_taxa
FROM uniprot
    INNER JOIN sequence_header ON
      sequence_header.uniprot_id = uniprot.id
    INNER JOIN tree_node uniprot_tn ON
      uniprot_tn.sequence_header_id = sequence_header.id
    INNER JOIN tree ON
      tree.id = uniprot_tn.tree_id
    INNER JOIN family ON
      family.canonical_tree_id = tree.id
    INNER JOIN tree_node all_tree_node ON
      all_tree_node.tree_id = tree.id
    INNER JOIN sequence_header all_sequence_header ON
      all_tree_node.sequence_header_id = all_sequence_header.id
    INNER JOIN uniprot all_uniprot ON
      all_sequence_header.uniprot_id = all_uniprot.id
WHERE
    uniprot.accession = %s AND
    family.family_type_id = 'G' AND
    family.status != 'bad' AND
    family.active = true
GROUP BY family.id
ORDER BY number_taxa DESC
LIMIT 1
)
"""
