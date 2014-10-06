BEGIN;
DROP TABLE metacyc.metacyc_version;
DROP TABLE metacyc.metacyc_node CASCADE;
DROP TABLE metacyc.metacyc_membership;

SET search_path=metacyc, public;

CREATE TABLE metacyc.metacyc_version (
	version text PRIMARY KEY
	);

CREATE TABLE metacyc.metacyc_node (
   id text PRIMARY KEY,
   common_name text,
   type text, 
   db text
);

CREATE TABLE metacyc.metacyc_protein ( 
	uniprot_accession text[]
) INHERITS (metacyc_node);

CREATE TABLE metacyc.metacyc_reaction (
	ec_number text,
	systematic_name text
) INHERITS (metacyc_node);

CREATE TABLE metacyc.metacyc_pathway ( ) INHERITS (metacyc_node);

CREATE TABLE metacyc.metacyc_gene ( ) INHERITS (metacyc_node);

CREATE TABLE metacyc.metacyc_membership (
   metacyc_member_id text NOT NULL,
   metacyc_group_id text NOT NULL,
   PRIMARY KEY (metacyc_member_id, metacyc_group_id)
);

\copy metacyc.metacyc_protein FROM '%(cur_dir)s/metacyc_protein.tab_delimited'
\copy metacyc.metacyc_gene FROM '%(cur_dir)s/metacyc_gene.tab_delimited'
\copy metacyc.metacyc_pathway FROM '%(cur_dir)s/metacyc_pathway.tab_delimited'
\copy metacyc.metacyc_reaction FROM '%(cur_dir)s/metacyc_reaction.tab_delimited'
\copy metacyc.metacyc_membership FROM '%(cur_dir)s/metacyc_membership.tab_delimited'
\copy metacyc.metacyc_version FROM '%(cur_dir)s/metacyc_version.tab_delimited'

COMMIT;
CREATE INDEX metacyc_prot_uniprot_accession on metacyc.metacyc_protein USING GIN (uniprot_accession);
CREATE INDEX metacyc_reaction_pkey on metacyc.metacyc_reaction (id);
CREATE INDEX metacyc_pathway_pkey on metacyc.metacyc_pathway (id);
CREATE INDEX metacyc_protein_pkey on metacyc.metacyc_protein (id);