/**
 * phylo_decoratetree - "Decorate" tree with links and
 * icons from the array of information passed as 
 * a JavaScript associative array.  
 * 
 * @author Christopher A. Meacham
 * Berkeley Phylogenomics Group
 * University of California, Berkeley
 * Copyright © 2008 Regents of the University of California
 *
 *
 */

var iconDir = '/static/img/icons/';
var showTaxon = true;
var showCommonName = false;
var showProtein = true;
var showUniProt = true;
var showUniProtId = false;
var showEnzymeCommission = true;
var showEvidence = true;
var showExperimentalOnly = true;
var showLiterature = true;
var showSubfamilyNodeNumber = false;
var showThreeDStructure = true;

var doColorSubfamilies = false;
var doSubfamilyCollapse = false;
var collapseThreshold = -1; // If greater than -1, will control setting of doSubfamilyCollapse based on number of leaves
var doColorOrtho = false;

var hideSubfamControls = false;
var hideOrthoControls = false;

var leftIdMap = null;
var subtreeLeftId = 1;
var highlightLeftIds = [];
var colorLeftIds = [];

function loadLeftIdMap(aRoot) {
	
	function loadLeftIdMap(aNode) {
		leftIdMap[aNode.leftId+''] = aNode;
//		aNode.text += "<br />"+aNode.leftId+" "+aNode.rightId;
		if (aNode.descendant) {
			loadLeftIdMap(aNode.descendant);
		}
		if (aNode.sibling) {
			loadLeftIdMap(aNode.sibling);
		}
	}

	leftIdMap = new Array();
	loadLeftIdMap(aRoot);
}

function isExperimentalCode(aCode) {
	return(aCode == 'EXP' || aCode == 'IDA' || aCode == 'IPI' || aCode == 'IMP' || 
		aCode == 'IGI' || aCode == 'IEP');
}

function initialCounts(aNode, aAnnotations, aCategory) {
// Counts are counts of sequences having evidence codes in a category.  These counts
//   are summed from leaves down for summarization at subfamily nodes. 
//   This function initializes the count for the category at this leaf node to one 
//   if the category is present.
	if (aAnnotations) {
		var _count = 0;
		for (i in aAnnotations) {
			aNode.evidenceCounts[aCategory] = 1;
			aNode.hasEvidence = 1;
			if (isExperimentalCode(aAnnotations[i][1])) {
				aNode.experimentalCounts[aCategory] = 1;
				aNode.hasExperimental = 1;	
				break;
			}
		}
	}
}
		
function collapseSubfamilies(aRoot, aBook) {
	
	function retrieveMRCA(aNode) {
		var queryStr = "";
		
		function retrieve(aNode) {
			if (aNode.taxid && aNode.taxid != "" && aNode.taxid != "N/A") {
				queryStr += "&taxid=" + aNode.taxid;
			}
			if (aNode.descendant) {
				retrieve(aNode.descendant);
			}
			if (aNode.sibling) {
				retrieve(aNode.sibling);
			}
		}
		
		retrieve(aNode.descendant); // Do not collect taxids of siblings of aNode
		queryStr = queryStr.substring(1);
		var ajaxObj = new sack();
		ajaxObj.node = aNode;
		ajaxObj.requestFile = '/rest/getMRCA?' + queryStr;
		ajaxObj.onCompletion = function() {
			var jsonObj = JSON.parse(ajaxObj.response);
			this.node.subfamilyObj = jsonObj;
			decorateSubfamilyNode(this.node, aBook);
		}
		ajaxObj.runAJAX();
	}
	
	function collapseCollapsableDescendantsOf(aNode) {
		var _n = aNode.descendant;
		while (_n) {
			if (_n.descendant && _n.subfamily_node != "") {
				_n.isCollapsed = true;
				_n.isSubfamilyNode = true;
				retrieveMRCA(_n); 
			}
			_n = _n.sibling;
		}
	}
	
	function addCounts(aCountArray1, aCountArray2) {
		if (aCountArray2) {
			for (key in aCountArray2) {
				if (!aCountArray1[key]) {
					aCountArray1[key] = 0;
				}
				aCountArray1[key] += aCountArray2[key];
			}
		}
	}
	 
	function compareSubfamilies(aNode) {
		if (aNode.descendant) {
			compareSubfamilies(aNode.descendant);
			aNode.subfamily_node = aNode.descendant.subfamily_node;
			aNode.subfamily_name = aNode.descendant.subfamily_name;
		    aNode.hasLiterature = aNode.descendant.hasLiterature;
		    aNode.hasEvidence = aNode.descendant.hasEvidence;
		    aNode.hasExperimental = aNode.descendant.hasExperimental;
		    aNode.hasSwissProt = aNode.descendant.hasSwissProt;
        aNode.hasThreeDStructure = aNode.descendant.hasThreeDStructure;
			aNode.evidenceCounts = [];
			addCounts(aNode.evidenceCounts, aNode.descendant.evidenceCounts);
			aNode.experimentalCounts = [];
			addCounts(aNode.experimentalCounts, aNode.descendant.experimentalCounts);
			if (aNode.subfamily_node && aNode.subfamily_node != "") {
				var n = aNode.descendant.sibling;
				while (n) {
					if (n.subfamily_node == aNode.subfamily_node) {
					    aNode.hasLiterature += n.hasLiterature;
              aNode.hasThreeDStructure += n.hasThreeDStructure;
					    aNode.hasEvidence += n.hasEvidence;
					    aNode.hasExperimental += n.hasExperimental;
					    aNode.hasSwissProt += n.hasSwissProt;
						addCounts(aNode.experimentalCounts, n.experimentalCounts);
						addCounts(aNode.evidenceCounts, n.evidenceCounts);
						n = n.sibling;
					} else {
						aNode.subfamily_node = "";
						collapseCollapsableDescendantsOf(aNode);
						break;
					}
				}
			
			} else {
				collapseCollapsableDescendantsOf(aNode);
			}
		}
		if (aNode.sibling) {
			compareSubfamilies(aNode.sibling);
		}
	}
		
	compareSubfamilies(aRoot);
}

function handleSubfamilyCollapse(doCollapse, aTree) {

	function expand(aRoot) {
		aRoot.tree.expandAll();
	}
	
	function collapse(aNode) {
		aNode.isCollapsed = (aNode.isSubfamilyNode == true); // aNode.isSubfamilyNode is either true or undefined
		if (aNode.descendant) {
			collapse(aNode.descendant);
		}
		if (aNode.sibling) {
			collapse(aNode.sibling);
		}
	}
	
	if (doCollapse && aTree && aTree.root) {
		collapse(aTree.root);
	} else {
		expand(aTree.root);
	}
}

function colorSubfamilySubtrees(aTree){
	var _count = 0;
	var _colors = ["red", "blue", "orange", "magenta", "cyan"];
	
	function findSubfamily(aNode) {
		if (aNode.isSubfamilyNode) {
			aTree.colorSubtree(aNode, _colors[_count % _colors.length]);
			_count++;
		} else {
			if (aNode.descendant) {
				findSubfamily(aNode.descendant);
			}
		}
		if (aNode.sibling) {
			findSubfamily(aNode.sibling);
		}
	}
	
	findSubfamily(aTree.root);
}

function handleColorSubfamilies(aDoColorSubfamilies, aTree) {
	if (aDoColorSubfamilies) {
		colorSubfamilySubtrees(aTree);
	} else {
		aTree.clearColorSubtrees();
	}
}

function colorOrtho(aTree, aOrthoObj, aLeftId) {
	if (aLeftId == 1) {
		var _baseTreeLeftId = aLeftId;
		var _baseTreeRightId = aTree.root.rightId + aLeftId - 1;
		//	alert(_baseTreeLeftId+" "+_baseTreeRightId);
		var _colors = ["seagreen", "purple", "tomato", "olive", "mediumvioletred"];
//	 	var _colors = ["red", "blue", "orange", "cyan", "magenta"];
		var _node = null;
		for (var i in aOrthoObj.superorthologous_left_ids) {
			_node = leftIdMap[aOrthoObj.superorthologous_left_ids[i]];
			//		if (_node.leftId >= _baseTreeLeftId && _node.leftId <= _baseTreeRightId) {
			aTree.colorSubtree(_node, _colors[i % _colors.length]);
		//		}
		}
	}
}

function colorLeftIdSubtrees(aTree, aLeftIdArray, aLeftId) {
	aTree.clearColorSubtrees();
	var _colorIndex = 0;
	var _color = "";
	var _baseTreeLeftId = aLeftId;
	var _baseTreeRightId = aTree.root.rightId + aLeftId - 1;
	var _colors = ["red", "blue", "orange", "cyan", "magenta"];
	var _node = null;
	for (var i in aLeftIdArray) {
		_node = leftIdMap[aLeftIdArray[i]];
		if (_node && _node.leftId >= _baseTreeLeftId && _node.leftId <= _baseTreeRightId) {
			_color = _colors[_colorIndex % _colors.length]
			if (_color == _node.lineDivStemColor) {
				_colorIndex++;
			}
			aTree.colorSubtree(_node, _color);
			_colorIndex++;
		}
	}
}

function handleColorOrtho(aDoColorOrtho, aTree) {
	if (aDoColorOrtho) {
		colorOrtho(aTree, orthoObj, subtreeLeftId);
	} else {
		aTree.clearColorSubtrees();
	}
}

function tip(aTipType) {
// Below, "this" refers to the containing HTML element (usually an <a>, <span>, or <img>), 
//   "parentNode" is the parent HTML <div> element that is the label in the tree,
//   and "node" is the node in the tree for which the <div> is the label.
	return(" onmouseover=\"showToolTip(this.parentNode.node, '" + aTipType + 
		"')\" onmouseout=\"UnTip()\" ");
}

function limitLength(aStr, aLength) {
	if (String(aStr).length > aLength) {
		aStr = String(aStr).substr(0, aLength) + "...";
	}
	return(aStr);
}

function validStr(aStr, aAlt) {
	if (aStr && aStr != "" && aStr != "N/A") {
		return(aStr);
	} else {
		return(aAlt)
	}
}

function decorateSubfamilyNode(aNode, aBook) {
	var _label = "";
	var _text = "";
	var _subfamilyText = "";
	if (aNode.subfamilyObj) {
		var commonName = aNode.subfamilyObj["common_name"];
		if (!commonName) {
			commonName = "[No&nbsp;common&nbsp;name]";
		}
		var scientificName = aNode.subfamilyObj["scientific_name"];
		if (!scientificName) {
			scientificName = "[No&nbsp;scientific&nbsp;name]";
		}
		_label += "<a href=\"http://www.uniprot.org/taxonomy/" + aNode.subfamilyObj["id"] + 
			"\"\n target=\"_blank\"" + tip("sf_common") + ">" + scientificName + "</a> ";
		_subfamilyText = limitLength(aNode.subfamily_name, 40);
		if (showSubfamilyNodeNumber) {
			_subfamilyText = " [" + aNode.subfamily_node + ": " + _subfamilyText + "]";
		} else {
			_subfamilyText = " [" + _subfamilyText + "]";
		}
		if (aBook && aBook != "") {
			_subfamilyText = "<a href=\"http://phylogenomics.berkeley.edu/book/book_info.php?book=" + aBook + 
				"&library=&subfamily=" + aNode.subfamily_node + "\" target=\"_blank\"" + 
				tip("sf_protein") + ">" + _subfamilyText + "</a>";
		}
		_label += _subfamilyText;
		if (showUniProt && aNode.hasSwissProt) {
			_label += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_swiss_flag_12.png\" onclick=\"this.parentNode.node.tree.collapseNode(this.parentNode.node)\"" + 
				tip("sf_swissprot") + "> ";
		}
		if ((showEvidence && aNode.hasEvidence && !showExperimentalOnly) || (showExperimentalOnly && aNode.hasExperimental)) {
			_label += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_flask_12.png\" onclick=\"this.parentNode.node.tree.collapseNode(this.parentNode.node)\"" + 
				tip("sf_evidence") + "> ";
		}
		if (showLiterature && aNode.hasLiterature) {
			_label += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_literature_12.png\" onclick=\"this.parentNode.node.tree.collapseNode(this.parentNode.node)\"" + 
				tip("sf_literature") + "> ";
		}
		if (showThreeDStructure && aNode.hasThreeDStructure) {
			_label += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_pdb_12.png\" onclick=\"this.parentNode.node.tree.collapseNode(this.parentNode.node)\"" + 
				tip("sf_structure") + "> ";
		}
		
		_text += scientificName;
		_text += "<br />" + commonName;
		_text += "<br />" + aNode.subfamily_node + ":&nbsp;"+ aNode.subfamily_name;
		if (aNode.hasSwissProt) {
			_text += "<br />SwissProt [" + aNode.hasSwissProt + "]";
		}
		if (aNode.hasEvidence) {
			_text += "<br />Annotation evidence [" + aNode.hasEvidence + "]";
		}
		if (aNode.hasExperimental) {
			_text += "<br />Experimental annotation evidence [" + aNode.hasExperimental + "]";
		}
		if (aNode.hasLiterature) {
			_text += "<br />Literature references [" + aNode.hasLiterature + "]";
		}
		if (aNode.hasThreeDStructure) {
			_text += "<br />PDB chain ids [" + aNode.hasThreeDStructure + "]";
		}
	}
	aNode.label = _label;
	aNode.text = _text;
	aNode.drawCollapsedLabelDiv(); // Force creation of aNode.collapsedDiv
}

function decorateNode(aNode, aTags, aBook) {
	if (aNode && aNode.subfamilyObj) {
		decorateSubfamilyNode(aNode, aBook);
	} else if (aNode && aTags) {
		aNode.hasSwissProt = 0;
		aNode.hasEvidence = 0;
		aNode.hasExperimental = 0;
		aNode.hasLiterature = 0;
    aNode.hasThreeDStructure = 0;
		var _htmlStr = "";
		var taxon = aTags["scientific_name"];
		var _protName = "";
		aNode.taxon = taxon;
		aNode.taxid = aTags["ncbi_taxid"];
		if (!aNode.descendant) {
			aNode.evidenceCounts = [];
			aNode.experimentalCounts = [];
			initialCounts(aNode, aTags.molecular_function, "function");
			initialCounts(aNode, aTags.biological_process, "process");
			initialCounts(aNode, aTags.cellular_component, "component");
		}
		if (showTaxon) {
			if (taxon != "N/A") {
				_htmlStr += "<a href=\"http://www.uniprot.org/taxonomy/" + aTags["ncbi_taxid"] + 
				"\"\n target=\"_blank\"" + tip("common") + ">" + limitLength(taxon, 20) + "</a> ";
			} else {
				_htmlStr += "Unknown taxon ";
			}
		}
		if (showCommonName && aTags["common_name"]) {
			var _cname = validStr(aTags["common_name"], "[No common name]");
			_htmlStr += "\n<a href=\"http://www.uniprot.org/taxonomy/" + aTags["ncbi_taxid"] + 
			"\"\n target=\"_blank\"" + tip("scientific") + ">" + _cname + "</a> ";
		}
		var protName = aTags["uniprot_de"];
		if (!protName) {
			protName = aTags["alignment_short_name"];
		}
		var _showingProtNameSection = (showProtein && protName && protName != "") || (showSubfamilyNodeNumber && aTags["subfamily_node"]);
		if (_showingProtNameSection) {
			_htmlStr += "\n [";
		}
		if (showSubfamilyNodeNumber && aTags["subfamily_node"]) {
			if (aBook && aBook != "") {
				_htmlStr += "<a href=\"http://phylogenomics.berkeley.edu/book/book_info.php?book=" + aBook + 
					"&library=gpcr&subfamily=" + aTags["subfamily_node"] + "\" target=\"_blank\"" + 
					tip("subfamily") + ">" + aTags["subfamily_node"] + "</a>";
			} else {
				_htmlStr += aTags["subfamily_node"];
			}
			if (showProtein && protName && protName != "") {
				_htmlStr += ": "; // Separate the two links
			}
		}
		if (showProtein && protName && protName != "") {
			_htmlStr += "<span "+ tip("protein") + ">" + limitLength(protName, 25) + "</span> ";
		}
		if (_showingProtNameSection) {
			_htmlStr += "]";
		}
		if (showUniProt && aTags["uniprot_accession"] != "") {
			if (aTags["uniprot_accession"] != "N/A") {
				_htmlStr += "\n<a href=\"http://www.uniprot.org/uniprot/" + 
					aTags["uniprot_accession"] + "\" target=\"_blank\"" + tip("uniprot") + ">";
			
				_htmlStr += aTags["uniprot_accession"];
				if (aTags["in_swissprot_f"] && aTags["in_swissprot_f"] == 1) {
					_htmlStr += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_swiss_flag_12.png\" alt=\"SwissProt\"> ";
				}
				_htmlStr += "</a> ";
			}
		}
		if (showUniProtId && aTags["uniprot_id"]) {
			_htmlStr += " " + aTags["uniprot_id"];
		}
		if (showEnzymeCommission && aTags["ec"] && aTags["ec"] != "") {
			_htmlStr += " "+aTags["ec"];
		}
		if ((showEvidence && aNode.hasEvidence && !showExperimentalOnly) || (showExperimentalOnly && aNode.hasExperimental)) {
			_htmlStr += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_flask_12.png\"" + tip("evidence") + "> ";
		}
		if (showLiterature && aTags["Literature"]) {
 			if (aTags["uniprot_accession"] && aTags["uniprot_accession"] != "N/A") {
				_htmlStr += "\n <a href=\"http://www.uniprot.org/uniprot/" + aTags["uniprot_accession"] + 
					"#section_ref\" target=\"_blank\"" + tip("literature") + ">" +
					" <img style=\"border: none;\" src=\"" + iconDir + "icon_literature_12.png\"></a> ";
			} else {
				_htmlStr += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_literature_12.png\"> ";
			}		
		}
    if (showThreeDStructure && aTags["ThreeDStructure"]) {
 			if (aTags["uniprot_accession"] && aTags["uniprot_accession"] != "N/A") {
				_htmlStr += "\n <a href=\"http://www.uniprot.org/uniprot/" + aTags["uniprot_accession"] + 
					"#section_x-ref_structure\" target=\"_blank\"" + tip("structure") 
          + ">" +
					" <img style=\"border: none;\" src=\"" + iconDir + "icon_pdb_12.png\"></a> ";
			} else {
				_htmlStr += "\n <img style=\"border: none;\" src=\"" + iconDir + "icon_pdb_12.png\"> ";
			}		
		}
		aNode.setLabelHTML(_htmlStr);
		
		var _textStr = "";
		if (aTags["in_swissprot_f"] && aTags["in_swissprot_f"] == 1) {
			aNode.hasSwissProt = 1;
		}
		if (aTags["scientific_name"]) {
			_textStr += aTags["scientific_name"];
		}
		if (aTags["common_name"]) {
			_textStr += "<br />" + aTags["common_name"];
		}
		if (protName) {
			_textStr += "<br />" + protName;
		}
		if (aTags["uniprot_accession"]) {
			_textStr += "<br />UniProt: " + aTags["uniprot_accession"];
		}
		if (aTags["uniprot_id"]) {
			_textStr += "<br />UniProt ID: " + aTags["uniprot_id"];
		}
		if (aTags["ec"] && aTags["ec"] != "") {
			_textStr += "<br />Enzyme Commission number: "+aTags["ec"];
		}
		if (aNode.hasEvidence) {
			_textStr += "<br />Annotation evidence";
		} 
		if (aNode.hasExperimental) {
			_textStr += "<br />Experimental annotation evidence";
		} 
		if (aTags["Literature"] && aTags["Literature"].length > 0) {
			aNode.hasLiterature = 1;
			_textStr += "<br />Literature available: " + aTags["Literature"].length + 
				" reference";
			if (aTags["Literature"].length > 1) {
				_textStr += "s";
			}
		}
    if (aTags["ThreeDStructure"] && aTags["ThreeDStructure"].length > 0) {
      aNode.hasThreeDStructure = 1;
      _textStr += "<br />ThreeD structure available: " + aTags["ThreeDStructure"].length
          + " PDB chain id";
      if (aTags["ThreeDStructure"].length > 1) {
        _textStr += "s";
      }
    }
		aNode.text = _textStr;
	}
	
}

function showToolTip(aNode, aTipType) {
	var _Str = "";
	var _goAbbreviations = [];
	var _annotationCategoryOrder = ["function", "process", "component"];
	
	function conCat(aStr1, aAlt1, aStr2, aAlt2, aSeparator) {
		aStr1 = validStr(aStr1, aAlt1);
		aStr2 = validStr(aStr2, aAlt2);
		if (aStr1 != "" && aStr2 != "") {
			return(aStr1 + aSeparator + aStr2);
		} else {
			return(aStr1 + aStr2);
		}
	}
	
	function showStr(aStr, aContinuationStr, aAlternateStr) {
		if (aStr && aStr != "" && aStr != "N/A") {
			phylo_Tip(aStr + aContinuationStr);
		} else {
			phylo_Tip(aAlternateStr);
		}
	}
	
	function show(aDataName, aContinuationStr, aAlternateStr) {
		showStr(jsonObj[aNode.label][aDataName], aContinuationStr, aAlternateStr);
	}
	
	function sf_show(aDataName, aContinuationStr, aAlternateStr) {
		showStr(aNode.subfamilyObj[aDataName], aContinuationStr, aAlternateStr);
	}
	
	function initializeGoAbbreviations() {
		_goAbbreviations.EXP = ["Inferred from Experiment"];
		_goAbbreviations.IDA = ["Inferred from Direct Assay"];
		_goAbbreviations.IPI = ["Inferred from Physical Interaction"];
		_goAbbreviations.IMP = ["Inferred from Mutant Phenotype"];
		_goAbbreviations.IGI = ["Inferred from Genetic Interaction"];
		_goAbbreviations.IEP = ["Inferred from Expression Pattern"];
		_goAbbreviations.ISS = ["Inferred from Sequence or Structural Similarity"];
		_goAbbreviations.ISO = ["Inferred from Sequence Orthology"];
		_goAbbreviations.ISA = ["Inferred from Sequence Alignment"];
		_goAbbreviations.ISM = ["Inferred from Sequence Model"];
		_goAbbreviations.IGC = ["Inferred from Genomic Context"];
		_goAbbreviations.RCA = ["Inferred from Reviewed Computational Analysis"];
		_goAbbreviations.TAS = ["Traceable Author Statement"];
		_goAbbreviations.NAS = ["Non-traceable Author Statement"];
		_goAbbreviations.IC = ["Inferred by Curator"];
		_goAbbreviations.ND = ["No biological Data available"];
		_goAbbreviations.IEA = ["Inferred from Electronic Annotation"];
		_goAbbreviations.NR = ["Not Recorded"];
		_goAbbreviations["N/A"] = ["Not Available"];
	}
	
	function concatenateGoAnnotations(aHeading, aAnnotations) {
		var _newStr = "";
		if (aAnnotations && aAnnotations.length > 0) {
			for (var i in aAnnotations) {
				if (isExperimentalCode(aAnnotations[i][1]) || !showExperimentalOnly) {
					_newStr += "&nbsp;&nbsp;&nbsp;&nbsp;" + aAnnotations[i][0] + " (" +
					aAnnotations[i][1] +
					")<br />";
					if (!_goAbbreviations[aAnnotations[i][1]]) {
						_goAbbreviations[aAnnotations[i][1]] = "Unknown abbreviation";
					}
					_goAbbreviations[aAnnotations[i][1]].present = true;
				}
			}
		}
		if (_newStr != "") {
			_Str += "&nbsp;" + aHeading + ":<br />" + _newStr;
		} 
	}
	
	function concatenateGoAbbreviations() {
		for (key in _goAbbreviations) {
			if (_goAbbreviations[key].present) {
				_Str += key + ": " + _goAbbreviations[key] + "<br />"
			}
		}
	}
	
	switch (aTipType) {
		case "common":
			_Str = conCat(jsonObj[aNode.label].scientific_name, "[No scientific name]", 
				jsonObj[aNode.label].common_name, "[No common name]", "<br />");
			showStr(_Str, "<br /><em>Click for UniProt taxonomy page...</em>", "[No taxonomic name]");
			break;
		case "scientific":
			_Str = conCat(jsonObj[aNode.label].scientific_name, "[No scientific name]", 
				jsonObj[aNode.label].common_name, "[No common name]", "<br />");
			showStr(_Str, "<br /><em>Click for UniProt taxonomy page...</em>", "[No taxonomic name]");
			break;
		case "subfamily":
			_Str = conCat(jsonObj[aNode.label].subfamily_node, "[No subfamily information]", 
				jsonObj[aNode.label].uniprot_de, "[No subfamily name]", ": ");
			showStr(_Str, "<br /><em>Click for PhyloFacts subfamily page...</em>", "[No subfamily information]");
			break;
		case "evidence":
			initializeGoAbbreviations();
			concatenateGoAnnotations("Molecular function", jsonObj[aNode.label].molecular_function);
			concatenateGoAnnotations("Biological process", jsonObj[aNode.label].biological_process);
			concatenateGoAnnotations("Cellular component", jsonObj[aNode.label].cellular_component);
			concatenateGoAbbreviations();
			if (_Str != "") {
				if (showExperimentalOnly) {
					_Str = "<em>Experimental annotation evidence</em><br />" + _Str;
				} else {
					_Str = "<em>Annotation evidence</em><br />" + _Str;
				}
			}
			showStr(_Str, "", "[No evidence information]");
			break;
		case "literature":
			_Str = jsonObj[aNode.label].Literature.length + " literature reference";
			if (jsonObj[aNode.label].Literature.length != 1) {
				_Str += "s";
			} 
			showStr(_Str, "<br /><em>Click for UniProt literature references...</em>", "[No literature information]")
			break;
    case "structure":
			_Str = jsonObj[aNode.label].ThreeDStructure.length + " PDB chain id";
			if (jsonObj[aNode.label].ThreeDStructure.length != 1) {
				_Str += "s";
			} 
			showStr(_Str, "<br /><em>Click for UniProt cross-references to the PDB...</em>", "[No ThreeD structure information]")
			break;
		case "protein":
			protName = jsonObj[aNode.label].uniprot_de;
			if (!protName) {
				protName = jsonObj[aNode.label].alignment_short_name;
			}
			if (jsonObj[aNode.label].subfamily_node) {
				protName = jsonObj[aNode.label].subfamily_node + ": " + protName;
			}
			showStr(protName, "", "[No available name]");
			break;
		case "uniprot":
			_Str = conCat(jsonObj[aNode.label].uniprot_accession, "[No UniProt information]", 
				jsonObj[aNode.label].uniprot_de, "[No UniProt descriptor]", ": ");
			showStr(_Str, "<br /><em>Click for UniProt page...</em>", "[No UniProt information]");
			break;
		case "sf_common":
			_Str = conCat(aNode.subfamilyObj.scientific_name, "[No scientific name]", 
				aNode.subfamilyObj.common_name, "[No common name]", "<br />");
			showStr(_Str, "<br /><em>Click for UniProt taxonomy page...</em>", "[No taxonomic name]");
			break;
		case "sf_swissprot":
			_Str = aNode.hasSwissProt + " SwissProt curated sequence";
			if (aNode.hasSwissProt != 1) {
				_Str += "s";
			}
			_Str += "<br /><em>Click to expand subfamily subtree...</em>";
			showStr(_Str, "", "[No SwissProt information]")
			break;
		case "sf_evidence":
			if (showExperimentalOnly) {
				if (aNode.experimentalCounts) {
					for (var i in _annotationCategoryOrder) {
						key = _annotationCategoryOrder[i];
					    if (aNode.experimentalCounts[key] && aNode.experimentalCounts[key] > 0) {
							_Str += "&nbsp;&nbsp;&nbsp;&nbsp;";
							switch (key) {
								case "function":
									_Str += "Molecular function";
									break;
								case "process":
									_Str += "Biological process";
									break;
								case "component":
									_Str += "Cellular component";
							}
							_Str += ": " + aNode.experimentalCounts[key] + "<br />";
						}
					}
				}
				if (_Str != "") {
					_Str = "<em>Experimental annotation evidence</em><br />&nbsp;Number of sequences with experimental evidence codes:<br />" + _Str;
				}
			} else {
				if (aNode.evidenceCounts) {
					for (var i in _annotationCategoryOrder) {
						key = _annotationCategoryOrder[i];
						if (aNode.evidenceCounts[key] && aNode.evidenceCounts[key] > 0) {
							_Str += "&nbsp;&nbsp;&nbsp;&nbsp;";
							switch (key) {
								case "function":
									_Str += "Molecular function";
									break;
								case "process":
									_Str += "Biological process";
									break;
								case "component":
									_Str += "Cellular component";
							}
							_Str += ": " + aNode.evidenceCounts[key] + "<br />";
						}
					}
				}
				if (_Str != "") {
					_Str = "<em>Annotation evidence</em><br />&nbsp;Number of sequences with evidence codes:<br />" + _Str;
				}
			}
			showStr(_Str, "<em>Click to expand subfamily subtree...</em>", "[No evidence information]")
			break;
		case "sf_literature":
			_Str = aNode.hasLiterature + " sequence";
			if (aNode.hasLiterature != 1) {
				_Str += "s have literature references";
			} else {
				_Str += " has a literature reference";
			}
			_Str += "<br /><em>Click to expand subfamily subtree...</em>";
			showStr(_Str, "", "[No literature information]")
			break;
    case "sf_structure":
			_Str = aNode.hasThreeDStructure + " sequence";
			if (aNode.hasThreeDStructure != 1) {
				_Str += "s have PDB chain ids";
			} else {
				_Str += " has a PDB chain id";
			}
			_Str += "<br /><em>Click to expand subfamily subtree...</em>";
			showStr(_Str, "", "[No ThreeD structure information]")
			break;
		case "sf_protein":
			showStr(aNode.subfamily_node + ": " + aNode.subfamily_name, "<br /><em>Click for PhyloFacts subfamily page...</em>", "[No available name]");
			break;
		default: phylo_Tip("'" + aTipType + "' is not available"); // Not normally called
	}
}

function phylo_Tip(aStr) {
	if (aStr && aStr != "") {
		Tip(aStr, FOLLOWMOUSE, false, WIDTH, -400, BORDERWIDTH, 2, PADDING, 8, DELAY, 250, FADEIN, 150, 
			FADEOUT, 0, DURATION, 0, FONTSIZE, '9pt');
	}
}


