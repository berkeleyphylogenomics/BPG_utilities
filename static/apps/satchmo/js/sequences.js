var theSeqDiv = null;
var nodes = [];
var allNodes = [];
var jsonSeqGlobal = null;
var jsonObjGlobal = null;
var displayString = new Array();
var leftCount = 0;

function addCharDiv(nodeIndex, stringIndex, aJsonSeq, id) {

	var color = "";

	switch(aJsonSeq[id][nodes[nodeIndex]].colors[stringIndex]){

		case "C":
			color = "cyan";
			break;
		case "B":
			color = "lightblue";
			break;
		case "G":
			color = "grey";
			break;
		case "N":
			color = "ghostwhite";
			break;
		default:
			color = "ghostwhite";
			break;
	}
	
	var topOffset = labelDivTop[nodes[nodeIndex]].replace("px", "");
	topOffset = parseInt(topOffset) + 340;
	topOffset = topOffset.toString() + 'px';
	//var topOffset = labelDivTop[nodes[nodeIndex]];
	var leftOffset = (document.body.clientWidth/4 + 4 + leftCount * 8) + 'px';
	//var leftOffset = (leftCount * 8) + 'px';
	leftCount++;
	//document.write(leftOffset);
	var retString = '<div style=\"font-size:11px;position:absolute;top:' + topOffset + ';left:' + leftOffset + ';background-color:' + color + ';font-family:\'Courier New\', Courier, monospace\">' +  aJsonSeq[id][nodes[nodeIndex]].aligned_sequence[stringIndex] + '</div>';

	//var retString = '<div style=\"font-size:11px;position:absolute;top:' + topOffset + ';background-color:' + color + ';font-family:\'Courier New\', Courier, monospace\">' +  aJsonSeq[id][nodes[nodeIndex]].aligned_sequence[stringIndex] + '</div>';

	return(retString);
}

function drawJsonSeq(jsonSeq, jsonObj, id) {
	var noOfNodes = -1;
	var i = 0;
	var nodeVal = null;
	var nodeNo = null;
	var tempStr = "";

	jsonSeqGlobal = jsonSeq;
	jsonObjGlobal = jsonObj;

	/*for(var key in labelDivTop){
		document.write("For " + key + ", we have " + labelDivTop[key]);
	}*/

	theSeqDiv = document.getElementById("seqdiv");

	noOfNodes = 0
	for(var key in jsonSeq[id]){
		nodes[noOfNodes] = key;
		//document.write(nodes[noOfNodes]);
		noOfNodes++;
	}
	noOfTotalNodes = newickTreeParseModified(jsonObj["__tree__"]);
	var temp = [];
	var count = 0;
	for(i = 0; i < noOfTotalNodes; i++){
		for(j = 0; j < noOfNodes; j++){
			if(allNodes[i] == nodes[j]){
				temp[count] = allNodes[i];
				count++;
				break;
			}
		}
	}

	for(i = 0; i < noOfNodes; i++){
		nodes[i] = temp[i];
	}

	for(i = 0; i < noOfNodes; i++){
		//document.write(nodes[i]);	
	}

	//tempStr += "<table width = \"100%\" border = \"5\">";
	//tempStr += "<tr>" + "<td>" + addCharDiv(0, 0, jsonSeq) + "</td><td>" + addCharDiv(0, 1, jsonSeq) + "</td></tr></table>";
	//theSeqDiv.innerHTML = "<table width = \"100%\" border = \"5\">" + "<tr>" + "<td>" + addCharDiv(0, 0, jsonSeq) + "</td><td>" + addCharDiv(0, 1, jsonSeq) + "</td></tr></table>";
	//theSeqDiv.innerHTML = tempStr;
	//tempStr += "<table width = \"100%\" border = \"0\" cellpadding = \"0\" cellspacing = \"0\">";
	//document.write("Hi" + jsonSeq[nodes[0]].colors[0]);
	//document.write(noOfNodes);
	for(i = 0; i < noOfNodes; i++){
		//tempStr += "<tr>";
		//jsonSeq[nodes[i]].aligned_sequence = jsonSeq[nodes[i]].aligned_sequence.replace('\n', '');
		//document.write(jsonSeq[nodes[i]].aligned_sequence.length);
		//document.write(jsonSeq[nodes[i]].aligned_sequence);
		for(j = 0; j < jsonSeq[id][nodes[i]].aligned_sequence.length; j++){
			//document.write("goin to addchardiv");
			//theSeqDiv.innerHTML += "<td>";
			//tempStr += "<td>" + addCharDiv(i, j, jsonSeq, id) + "</td>";
			tempStr += addCharDiv(i, j, jsonSeq, id);
 		}
		//theSeqDiv.innerHTML += jsonSeq[nodes[i]].aligned_sequence + "<br />";
		//theSeqDiv.innerHTML += "<br />";
		//addCharDiv(0, 0, jsonSeq);
		//document.write(theSeqDiv.innerHTML);
		//tempStr += "</tr>";
		leftCount = 0;
	}
	//tempStr += "</table>";
	displayString[id] = tempStr;
	//theSeqDiv.innerHTML = tempStr;
}

function drawAllJsonSeq(jsonSeq, jsonObj) {

	for(var key in jsonSeq){
		drawJsonSeq(jsonSeq, jsonObj, key);
	}
	theSeqDiv.innerHTML = displayString[1];
}

function newickTreeParseModified(treeStr) {
	var root = null;
	var ch = " ";
	var chIndex = -1;
	var counter = 0;

	function getCh() {
		do {
			chIndex++;	
		} while (chIndex<treeStr.length && treeStr.charAt(chIndex)=="\n" &&
			treeStr.charAt(chIndex)=="\r" && treeStr.charAt(chIndex)=="\t");
		if (chIndex<treeStr.length) {
			ch = treeStr.charAt(chIndex);
		} else {
			ch = ")";  // Close unmatched open parentheses
		}
	}
	
	function processNodeLabel(aNode) {
		var labelStr = "";
		do {
			labelStr += ch;
			getCh();
		} while (!(ch==" " || ch=="(" || ch==")" || ch=="," || ch==":" || ch=="[" || ch==";"));
		//labelStr = labelStr.replace(/_/g, " ");
		if (labelStr.indexOf("*LOST") > -1) {
			aNode.label = labelStr;
		}
		else {
			aNode.label = labelStr;
			aNode.text = labelStr;
			allNodes[counter] = labelStr;
			counter++;
		}
	}
	
	function processNodeLength(aNode) {
		var lengthStr = "";
		do {
			getCh();
		} while (ch==" ");
		while ((ch>="0" && ch<="9") || ch=="." || ch=="+" || ch=="-" || ch=="E" || ch=="e") {
			lengthStr += ch;
			getCh();
		}
		aNode.length = (+lengthStr);  // Force conversion of string to number		
	}
	
	function processBracketStr(aStr, aNode) {
		var position = aStr.indexOf("&&NHX");
		if (position>-1) {
			nhx_isPresent = true;
			aStr = aStr.substring(position+5);
			if (aStr.charAt(0)==":") {
				aStr = aStr.substring(1);
				var phrase = aStr.split(":");
				for (var i=0; i<phrase.length; i++) {
					var term = phrase[i].split("=");
					switch (term[0]) {
						case "S":
							aNode.nhx_S = term[1];
							break;
						case "D":
							aNode.nhx_D = term[1];
							break;
						case "B":
							aNode.nhx_B = term[1];
							aNode.length = (+aNode.nhx_B);
							break;
						case "Co":
							if (term[1]="Y") {
								aNode.isCollapsed = true;
							}
					}
				}
			}
		} else {
			aNode.label = aStr;
		}
	}
	
	function processBetweenBrackets(aNode) {
		var bracketStr = "";
		var ofStr = "";
		do {
			getCh();
		} while (ch==" ");
		if (ch >= "0" && ch <= "9") {
			while (ch != "]") {
				bracketStr += ch;
				getCh();
			}
			getCh();
			while (ch != "[") {
				ofStr += ch;
				getCh();
			}
			aNode.label = "["+bracketStr+"]"+ofStr;
			ch = " ";
			chIndex--;  // Back up so that "[" will be read again by ProcessNode()
//alert("Lost: "+bracketStr);			
		} else {
			while (ch != "]") {
				bracketStr += ch;
				getCh();
			}
			getCh();
			processBracketStr(bracketStr, aNode);
//alert("Bracket: "+bracketStr);			
		}
	}
	
	function processNode(ancNode) {
		var doingDescendants = false;
		var node = new newickNode();
		node.ancestor = ancNode;
		if (node.ancestor) {
			node.leftId = node.ancestor.rightId;
			node.rightId = node.leftId+1;
			node.ancestor.rightId = node.rightId+1;
		} else {
			node.leftId = 1;
			node.rightId = 2;
		}
		while (ch==" ") {
			getCh();
		}
		while (ch!=";" && ((ch!="," && ch!=")") || doingDescendants)) {
			switch (ch) {
				case "(":
					doingDescendants = true;
					getCh();
					node.descendant = processNode(node);
					node.lastDescendant = node.descendant;
					node.rightId = node.lastDescendant.rightId+1;
					break;
				case ",":
					getCh();
					node.lastDescendant.sibling = processNode(node);
					node.lastDescendant = node.lastDescendant.sibling;
					node.rightId = node.lastDescendant.rightId+1;
					break;
				case ")":
					doingDescendants = false;
					do {
						getCh();
					} while (ch==" ");
					break;
				case ":":
					processNodeLength(node);
					break;
				case "[":
					processBetweenBrackets(node);
					break;
				default:
					if (!(ch==";" || ch=="," || ch=="(" || ch==")")) {
						processNodeLabel(node);
					}
			}
			while (ch==" ") {
				getCh();
			}
			
		}
		return (node);
	}

	function labelNode(aNode) {
		aNode.label = aNode.leftId+"  "+aNode.rightId;
		if (aNode.descendant) {
			labelNode(aNode.descendant);
		}
		if (aNode.sibling) {
			labelNode(aNode.sibling);
		}
	}
	
	root = processNode(null);
	if (root) {
//		labelNode(root);
	}
	return(counter);
}
