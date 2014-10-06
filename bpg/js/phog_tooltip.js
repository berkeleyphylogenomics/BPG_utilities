/**
 *phog_tooltip.js - Tooltip display for PHOG pages.
 * 
 * @author Christopher A. Meacham
 */

var inputTip = 'Accepted inputs &mdash; <strong>for protein sequences only!</strong> &mdash; include:'+
	'<ul><li>UniProt accession (e.g., P28222)</li>'+
	'<li>UniProt ID (e.g., 5HT1B_HUMAN)</li>'+
	'<li>GenBank ID (e.g., 112821)</li>'+
	'<li>PHOG accession (e.g., PHOG006510081)</li></ul>'+
	'Disallowed inputs: Please note that we do not accept GenBank accessions '+
	'(e.g., CAC82475),  or identifiers from other databases (e.g., PDB). We also '+
	'do not accept IDs or accessions corresponding to nucleotide sequences.<br /><br />'+
	'GenBank accessions and IDs can be differentiated as follows: GenBank '+
	'accessions are alpha-numeric whereas GenBank IDs are strictly numeric.';	

/* This one excludes PHOG accessions */
var netscope_inputTip = 'Accepted inputs &mdash; <strong>for protein sequences only!</strong> &mdash; include:'+
	'<ul><li>UniProt accession (e.g., P28222)</li>'+
	'<li>UniProt ID (e.g., 5HT1B_HUMAN)</li>'+
	'<li>GenBank ID (e.g., 112821)</li>'+
	'Disallowed inputs: Please note that we do not accept GenBank accessions '+
	'(e.g., CAC82475),  or identifiers from other databases (e.g., PDB). We also '+
	'do not accept IDs or accessions corresponding to nucleotide sequences.<br /><br />'+
	'GenBank accessions and IDs can be differentiated as follows: GenBank '+
	'accessions are alpha-numeric whereas GenBank IDs are strictly numeric.';	

function showInputTip() {
	bpgToolTip(inputTip);
}

function showNetscopeInputTip() {
	bpgToolTip(netscope_inputTip);
}

function showNSeqsTip(aNSeqs) {
	bpgToolTip("<em>Click to view the full family tree (" + aNSeqs + " sequences)...</em>");
}

/**
*   The convention followed here is that for a datagrid cell the tip type is the (internal) column
*     name and for a datagrid column heading the tip type is the (internal) column name preceeded 
*	  by an underscore character.  For example, the column labelled "UniProt ID", which is named
*     "uniprot_id", has a tip type of "_uniprot_id" for the text of the column heading and 
*     "uniprot_id" for the tip type of the entries (cells) in the datagrid column.
*  
*/

function showToolTip(aTipType) {
	
	function showStr(aStr) {
		bpgToolTip(aStr);
	}
	
	switch (aTipType) {
		
		case "_aligned_length":
			showStr("Number of residues in the BLAST aligned region");
			break;
		case "_biocyc":
			showStr("Links to the BioCyc collection of metabolic pathways databases");
			break;
		case "biocyc":
			showStr("<em>Click for BioCyc metabolic pathway page...</em>");
			break;
		case "_blast_alignment":
			showStr("Links to pairwise BLAST alignments of subject sequence with query sequence");
			break;
		case "blast_alignment":
			showStr("<em>Click to view pairwise BLAST alignment of this subject sequence with query sequence...</em>");
			break;
		case "_blast_score":
			showStr("BLAST score of the highest-scoring aligned region for each sequence");
			break;
		case "_book":
			showStr("Links to PhyloFacts protein family books");
			break;
		case "book":
			showStr("<em>Click for PhyloFacts protein family book...</em>");
			break;
		case "_book_type":
			showStr("Alignment type &ndash; global or local");
			break;
		case "book_type":
			showStr("Alignment type");
			break;
		case "_description":
			showStr("Protein description");
			break;
		case "description":
			showStr("Protein description");
			break;
		case "dip":
			showStr("<em>Links to the Database of Interacting Proteins </em>");
			break;
		case "_dip":
			showStr("<em>Click for DIP source page...</em>");
			break;
		case "_e_value":
			showStr("BLAST E-value &ndash; the number of hits like this one expected at random in a database of this size");
			break;
		case "_ec":
			showStr("Links to the Enzyme Commission classification");
			break;
		case "ec":
			showStr("<em>Click for Enzyme Commission page...</em>");
			break;
		case "_experimental_evidence":
			showStr("Links to Gene Ontology experimental evidence for biological process, molecular function and cellular component");
			break;
		case "experimental_evidence":
			showStr("<em>Click for Gene Onotology evidence page...</em>");
			break;
		case "full_tree":
			showStr("<em>Click to view the full family tree...</em>");
			break;
		case "_gene":  
			showStr("Gene ID/accession linked to the UniProt page");
			break;
		case "gene":   
			showStr("<em>Click for UniProt page...</em>");
			break;
		case "_go_evidence":  
			showStr("Gene Ontology evidence for biological process, molecular function and cellular component");
			break;
		case "go_evidence":
			showStr("<em>Click for UniProt GO evidence page...</em>");
			break;
		case "_in_swissprot":
			showStr("Links for SwissProt curated sequences");
			break;
		case "in_swissprot":
			showStr("<em>Click for SwissProt entry...</em>");
			break;
		case "_kegg":
			showStr("Links to the Kyoto Encyclopedia of Genes and Genomes");
			break;
		case "kegg":
			showStr("<em>Click for KEGG page...</em>");
			break;
		case "_literature":
			showStr("Links to UniProt literature references");
			break;
		case "literature":
			showStr("<em>Click for UniProt literature references...</em>");
			break;
		case "ncbi":   
			showStr("<em>Click for NCBI page...</em>");
			break;
	        case "netscope":
		        showStr("Following this link will open our PHOG PPI resource, with a graphical display of<ul><li>observed interacting partners (currently taken from UCLA's Database of Interacting Proteins), and</li><br ><li>predicted interacting partners, found by interolog analysis using the PHOG algorithm.</li></ul>");
			break;
		case "_num_seqs":
			showStr("Number of non-redundant sequences");
			break;
		case "_orthologs":
			showStr("Links to the PhyloFacts Orthology Group pages listing orthologs of these sequences");
			break;
		case "orthologs":
			showStr("<em>Click to view the PhyloFacts Orthology Group page listing orthologs for this sequence...</em>");
			break;
		case "_percent_cover":
			showStr("Percent coverage");
			break;
		case "_percent_id":
			showStr("Percent identity of amino acid residues with query sequence over the aligned region");
			break;
		case "_pfam":
			showStr("Pfam designation of domains found in sequence linked to the Pfam page");
			break;
		case "pfam":
			showStr("<em>Click for Pfam page...</em>");
			break;
		case "_phog":
			showStr("Links to PhyloFacts Orthology Group pages");
			break;
		case "phog":
			showStr("<em>Click for PHOG page...</em>");
			break;
		case "_phog_threshold":
			showStr("Larger thresholds increase sensitivity, but can take more time.");
			break;
		case "_phog_tree":
			showStr("Links to PhyloFacts trees for the full family and for the PhyloFacts Orthology Group");
			break;
		case "phog_tree":
			showStr("<em>Click to view PHOG tree...</em>");
			break;
		case "phyloscope":
			showStr("<em>Click to view tree with PhyloScope...</em>");
			break;
		case "_ppi":
			showStr("Links to information on Protein-Protein Interactions");
			break;
		case "ppi":
			showStr("<em>Click for ...</em>");
			break;
		case "_species":
			showStr("Scientific (and common) name linked to UniProt taxonomy page");
			break;
		case "species":
			showStr("<em>Click for UniProt taxonomy page...</em>");
			break;
		case "_uniprot_id":
			showStr("UniProt ID");
			break;
		case "uniprot_id":
			showStr("UniProt ID");
			break;
		case "uniprot_id":
			showStr("UniProt accession");
			break;
		default: bpgToolTip("'" + aTipType + "' is not a defined tool tip type"); // For debugging; not normally called
	}
}
