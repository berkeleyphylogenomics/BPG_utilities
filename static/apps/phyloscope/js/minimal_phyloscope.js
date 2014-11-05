var done = false;
var jsonObj = null;
var tree = null;
var theTreeDiv = null;

function message(aMessage) {
	theMessageDiv.innerHTML = aMessage;
}

function doRedraw() {
  if (tree) {
    tree.draw();
  }
}

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

function colorLeftIdSubtrees(aTree, aLeftIdArray, aLeftId) {
	aTree.clearColorSubtrees();
	var _colorIndex = 0;
	var _color = "";
	var _baseTreeLeftId = aLeftId;
	var _node = null;
  _node = leftIdMap[aLeftId];
  var _baseTreeRightId = _node.rightId;
	var _colors = ["red", "blue", "orange", "cyan", "magenta"];
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
	subtreeLeftId = aSubtreeLeftId;
	highlightLeftIds = aHighlightLeftIds;
	theTreeDiv = document.getElementById("treediv"); 
	
	theTreeDiv.innerHTML = "<br /><br /><br />&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;Loading tree...";
}

function drawTreeStr(aTreeStr) {
	theTreeDiv.innerHTML = "<br /><br /><br />&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;Loading tree...";
	var root = newickTreeParse(aTreeStr);

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

