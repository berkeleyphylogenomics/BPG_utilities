/**
 * phyloscope - Contains most of the top-level functions that interact with 
 * the PhyloScope page DOM objects.  
 *
 * @author Christopher A. Meacham
 * Berkeley Phylogenomics Group
 * University of California, Berkeley
 * Copyright © 2008 Regents of the University of California
 *
 *
 */

var familyAccessionNumber = "";
var checkBoxesChanging = false;
var done = false;
var jsonObj = null;
var orthoObj = null;
var tree = null;
var theTreeDiv = null;

function message(aMessage) {
	theMessageDiv.innerHTML = aMessage;
}

// Used to redraw tree when display options have changed
function doRedraw() {
  if (tree) {
    traverseAll(tree.root, loadTags);
    tree.draw();
  }
}


// Sets the checkboxes based on the boolean control values in phylo_decoratetree.js
function setCheckBoxes() {
	checkBoxesChanging = true;
  	document.getElementById('cbTaxon').checked = showTaxon;
  	document.getElementById('cbCommonName').checked = showCommonName;
  	document.getElementById('cbProtein').checked = showProtein;
  	document.getElementById('cbUniProt').checked = showUniProt;
  	document.getElementById('cbUniProtId').checked = showUniProtId;
  	document.getElementById('cbEnzymeCommission').checked = showEnzymeCommission;
  	document.getElementById('cbEvidence').checked = showEvidence;
	if (!showEvidence) {
		showExperimentalOnly = false;
	}
  	document.getElementById('cbExperimentalOnly').checked = showExperimentalOnly;
  	document.getElementById('cbExperimentalOnly').disabled = !showEvidence;
	if (document.getElementById('cbExperimentalOnly').disabled) {
		document.getElementById('cbExperimentalOnlyText').style.color = "gray";
	} else {
		document.getElementById('cbExperimentalOnlyText').style.color = "black";
	}
  	document.getElementById('cbLiterature').checked = showLiterature;
  	document.getElementById('cbColorOrtho').checked = doColorOrtho;
	checkBoxesChanging = false;
}


// Set the boolean control variables in phylo_decoratetree.js from the state of the check boxes
function cbChanged() {
  if (!checkBoxesChanging) {
  	showTaxon = document.getElementById('cbTaxon').checked;
  	showCommonName = document.getElementById('cbCommonName').checked;
  	showProtein = document.getElementById('cbProtein').checked;
  	showUniProt = document.getElementById('cbUniProt').checked;
  	showUniProtId = document.getElementById('cbUniProtId').checked;
  	showEnzymeCommission = document.getElementById('cbEnzymeCommission').checked;
  	showEvidence = document.getElementById('cbEvidence').checked;
	if (!showEvidence) {
		document.getElementById('cbExperimentalOnly').checked = false;
	}
  	showExperimentalOnly = document.getElementById('cbExperimentalOnly').checked;
  	document.getElementById('cbExperimentalOnly').disabled = !showEvidence;
	if (document.getElementById('cbExperimentalOnly').disabled) {
		document.getElementById('cbExperimentalOnlyText').style.color = "gray";
	} else {
		document.getElementById('cbExperimentalOnlyText').style.color = "black";
	}
  	showLiterature = document.getElementById('cbLiterature').checked;
  	tree.loadLabelDivs();
  	doRedraw();
  }
}

// Handle toggling of coloring of super-orthologous groups
function cbColorOrthoChanged() {
  if (!checkBoxesChanging && tree) {
  	doColorOrtho = document.getElementById('cbColorOrtho').checked;
	handleColorOrtho(doColorOrtho, tree);
  }
}

// Generate the HTML to display the node information based on the control 
//  variables in phylo_decoratetree.js
function loadTags(aNode) {
	var treeTag = aNode.label;
	if (jsonObj[treeTag]) {
		decorateNode(aNode, jsonObj[treeTag], familyAccessionNumber);
	} else {
//		aNode.labelDiv.innerHTML = "Tree tag not found: " + treeTag;
//		aNode.labelDiv.style.color = "red";
	}
}

function htmlEntityDecode(str) {
	try {
		var tarea = document.createElement('textarea');
		tarea.innerHTML = str; 
		return tarea.value;
		tarea.parentNode.removeChild(tarea);
	} catch(e) {
//for IE add <div id="htmlconverter" style="display:none;"></div> to the page
		document.getElementById("htmlconverter").innerHTML = '<textarea id="innerConverter">' + str + '</textarea>';
		var content = document.getElementById("innerConverter").value;
		document.getElementById("htmlconverter").innerHTML = "";
		return content;
	}
}

// This is the old draw method, called from the AJAX object; no longer used
function draw(){
  theTreeDiv.innerHTML = "";
  var root = newickTreeParse(jsonObj["__tree__"]);
  
  if (subtreeLeftId != 1) {
  	treeStr = newickTreeWriteSubtree(newickTreeLeftIdNode(root, subtreeLeftId));
 	root = newickTreeParse(treeStr);
	if (root) {
		root.length = 0.0;
		newickTreeAdjustNodeIds(root, subtreeLeftId);
	}
  }
    
  tree = new newickTreeGraphic(root, 'treediv');
  loadLeftIdMap(tree.root);
  tree.setPercentOffsetAndWidth(0, 100);
  tree.onNodeButtonMouseOver = function(aNode) {
//    phylo_Tip(aNode.text);
  }
  tree.onNodeButtonMouseOut = function(aNode) {
    UnTip();
  }
  tree.onNodeMouseOver = function(aNode) {
  	// Do nothing - necessary to override display of aNode.text
  }
  traverseLabelled(tree.root, loadTags);
  tree.draw();
  for (var i in highlightLeftIds) {
  	var node = leftIdMap[highlightLeftIds[i]];
	if (node && !node.descendant) {
	  node.isInSet = true;
	}
  }
  tree.drawHighlightNodeSet();
}

function initControls() {
	if (hideOrthoControls) {
  		document.getElementById('divColorOrtho').style.display = "none";
	}	
}

function drawJsonObj(aJsonObj) {
	theTreeDiv.innerHTML = "<br /><br /><br />&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;Loading tree...";
	var root = newickTreeParse(aJsonObj["__tree__"]);

	if (subtreeLeftId != 1) {
		var treeStr = newickTreeWriteSubtree(newickTreeLeftIdNode(root, subtreeLeftId));
		root = newickTreeParse(treeStr);
		if (root) {
			root.length = 0.0;
			newickTreeAdjustNodeIds(root, subtreeLeftId);
		}
	}

	tree = new newickTreeGraphic(root, 'treediv');
	loadLeftIdMap(tree.root);

	tree.setPercentOffsetAndWidth(0, 100);
	tree.onNodeMouseOver = function(aNode) {
		// Do nothing - necessary to override display of aNode.text
	}
	tree.onNodeButtonMouseOver = function(aNode) {
		phylo_Tip(aNode.text);
	}
	tree.onNodeButtonMouseOut = function(aNode) {
		UnTip();
	}
  
	theTreeDiv.innerHTML = "";
	traverseLabelled(tree.root, loadTags);

	if (orthoObj != null) {
		colorOrtho(tree, orthoObj, subtreeLeftId);
	}

	if (collapseThreshold > -1) {
		countTips(tree.root);
	}

	setCheckBoxes();

	tree.draw();
	
	if (highlightLeftIds.length > 0 || colorLeftIds.length > 0) {
		loadLeftIdMap(tree.root);
	}

	if (colorLeftIds.length > 0) {
		colorLeftIdSubtrees(tree, colorLeftIds, subtreeLeftId);
	}

	if (highlightLeftIds.length > 0) {
		for (var i in highlightLeftIds) {
			var node = leftIdMap[highlightLeftIds[i]];
			if (node && !node.descendant) {
				node.isInSet = true;
			}
		}
		tree.drawHighlightNodeSet();
	}

}

