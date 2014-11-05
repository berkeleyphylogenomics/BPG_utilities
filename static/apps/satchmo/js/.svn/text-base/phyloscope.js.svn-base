var bookAccessionNumber = "";
var checkBoxesChanging = false;
var done = false;
var jsonObj = null;
var orthoObj = null;
var tree = null;
var theTreeDiv = null;

function message(aMessage) {
	theMessageDiv.innerHTML = aMessage;
}

function doRedraw() {
  if (tree) {
    traverseAll(tree.root, loadTags);
    tree.draw();
  }
}

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
  	document.getElementById('cbSubfamilyNodeNumber').checked = showSubfamilyNodeNumber;
  	document.getElementById('cbColorSubfamilies').checked = doColorSubfamilies;
  	document.getElementById('cbColorOrtho').checked = doColorOrtho;
  	document.getElementById('cbCollapseSubfamilies').checked = doSubfamilyCollapse;
	checkBoxesChanging = false;
}

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
  	showSubfamilyNodeNumber = document.getElementById('cbSubfamilyNodeNumber').checked;
  	tree.loadLabelDivs();
  	doRedraw();
  }
}

function cbCollapseSubfamiliesChanged () {
  if (!checkBoxesChanging && tree) {
  	doSubfamilyCollapse = document.getElementById('cbCollapseSubfamilies').checked;
	handleSubfamilyCollapse(doSubfamilyCollapse, tree);
	doRedraw();
  }
}

function cbColorSubfamiliesChanged () {
  if (!checkBoxesChanging && tree) {
  	doColorSubfamilies = document.getElementById('cbColorSubfamilies').checked;
	handleColorSubfamilies(doColorSubfamilies, tree);
  }
}

function cbColorOrthoChanged () {
  if (!checkBoxesChanging && tree) {
  	doColorOrtho = document.getElementById('cbColorOrtho').checked;
	handleColorOrtho(doColorOrtho, tree);
  }
}

function loadTags(aNode) {
	var treeTag = aNode.label;
	if (jsonObj[treeTag] || aNode.subfamilyObj) {
		decorateNode(aNode, jsonObj[treeTag], bookAccessionNumber);
	} else {
//		aNode.labelDiv.innerHTML = "Tree tag not found: " + treeTag;
//		aNode.labelDiv.style.color = "red";
	}
}

function loadSubfamily(aNode) {
	var treeTag = aNode.label;
	if (jsonObj[treeTag]) {
		if (jsonObj[treeTag].subfamily_node) {
			aNode.subfamily_node = jsonObj[treeTag].subfamily_node;
		} else {
			aNode.subfamily_node = "";
		}
		if (jsonObj[treeTag].alignment_uniprot_de_name && jsonObj[treeTag].alignment_uniprot_de_name != "N/A") {
			aNode.subfamily_name = jsonObj[treeTag].alignment_uniprot_de_name;
		} else if (jsonObj[treeTag].alignment_short_name && jsonObj[treeTag].alignment_short_name != "N/A") {
			aNode.subfamily_name = jsonObj[treeTag].alignment_short_name;
		} else if (jsonObj[treeTag].uniprot_de && jsonObj[treeTag].uniprot_de != "N/A") {
			aNode.subfamily_name = jsonObj[treeTag].uniprot_de;
		} else {
			aNode.subfamily_name = "[No subfamily name]";
		}
	}
}

//Called whenever a JSON object is to be interpreted in a template
//Make a text area, puts the JSON string as its value, returns that 
//value and finally removes the text area.
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
  traverseLabelled(tree.root, loadSubfamily);
//  collapseSubfamilies(tree.root, bookAccessionNumber);
  if (doColorSubfamilies) {
  	colorSubfamilySubtrees(tree)
  }
//  handleSubfamilyCollapse(doSubfamilyCollapse, tree);
  tree.draw();
  for (var i in highlightLeftIds) {
  	var node = leftIdMap[highlightLeftIds[i]];
	if (node && !node.descendant) {
	  node.isInSet = true;
	}
  }
  tree.drawHighlightNodeSet();
}

function initializePhyloScope(aBook, aSubtreeLeftId, aHighlightLeftIds) {
//	if (aHighlightLeftIds) {alert(aHighlightLeftIds.join("-"))}
	bookAccessionNumber = aBook;
	subtreeLeftId = aSubtreeLeftId;
	highlightLeftIds = aHighlightLeftIds;
	theTreeDiv = document.getElementById("treediv"); 
	var ajaxURL = 'http://phylogenomics.berkeley.edu/rest/index/annotateTree?accession=' + aBook;
	ajaxObj = new sack(); 
	ajaxObj.requestFile = ajaxURL;
	ajaxObj.onCompletion = function() {
		var jsonStr = ajaxObj.response;
		jsonObj = JSON.parse(jsonStr);
		draw(jsonObj);

		var ajaxURL_ortho = 'http://phylogenomics.berkeley.edu/rest/index/getSuperorthologousNodes?accession=' + aBook;
		ajaxObj2 = new sack();
		ajaxObj2.requestFile = ajaxURL_ortho;
		ajaxObj2.onCompletion = function() {
			var jsonStr2 = ajaxObj2.response;
			orthoObj = JSON.parse(jsonStr2);
			colorOrtho(tree, orthoObj, subtreeLeftId);			
			done = true;
		}
		ajaxObj2.runAJAX();
	}
	
	theTreeDiv.innerHTML = "<br /><br /><br />&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;Loading tree...";
	setCheckBoxes();
	ajaxObj.runAJAX();
}

function drawJsonObj(aJsonObj) {
	theTreeDiv.innerHTML = "<br /><br /><br />&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;Loading tree...";
	var root = newickTreeParse(aJsonObj["__tree__"]);

	if (subtreeLeftId != 1) {
		doSubfamilyCollapse = false;
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
	traverseLabelled(tree.root, loadSubfamily);

	if (orthoObj != null) {
		colorOrtho(tree, orthoObj, subtreeLeftId);
		doColorSubfamilies = false;
		doCollapseSubfamiles = false;
	}

	collapseSubfamilies(tree.root, bookAccessionNumber);
	if (doColorSubfamilies) {
		colorSubfamilySubtrees(tree)
	}
	handleSubfamilyCollapse(doSubfamilyCollapse, tree);

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

