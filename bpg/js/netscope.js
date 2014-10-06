
var unique_phog_list = new Array();
var unique_phog_count = 1;
var edges;
var edge_count = 0;
//var node_list = new Array();
var node_count = 0;
var query = "";
var taxon_scientific_name = "";
var taxon_common_name = "";
var ortho_type = "";

function nodeClick(node) {
  //alert("Node click not implemented.  This will open a sequence-specific BPG page.");
  //window.open("http://phylofacts.berkeley.edu/seq/"+node+"/");
}

function nodeTip(node,description) {
  var tip = "<strong>Protein description:</strong>,<br /> &nbsp; &nbsp; &middot; <a href='http://www.uniprot.org/uniprot/"+node+"' target='_blank' title='Click for UniProt page...'>"+description+"</a><br />";

  tip += "<br />Click <a href='/phog/net/?sequence_id="+node+"&ortholog_type="+ortho_type;

  tip += "'>here</a> to re-center on this node...";
  tip_title = ""+node;
  Tip(tip, STICKY, true, CLOSEBTN, true, FONTSIZE, '9pt', OFFSETX, 4, OFFSETY, 4, TITLE, tip_title);

}

function PHOGNodeTip(node) {
  var genes = node_list[node]['genes_from_taxon'];
  var description = node_list[node]['description'];
  var tip = "<strong>Protein description:</strong>,<br /> &nbsp; &nbsp; &middot; <a href='http://phylofacts.berkeley.edu/phog/"+node+"/' target='_blank' title='Click for orthologs...'>"+node_list[node]['description']+"</a>,<br />";

  if (genes) {
    if (genes.length > 0) {
      tip+=",<br /><strong>Genes from query taxon:</strong>,<br />";
      for (var i = 0; i < genes.length; i++) {
	tip += " &nbsp; &nbsp; &middot; <a href='http://www.uniprot.org/uniprot/"+genes[i];
	tip += "' target='_blank' title='Click for UniProt page...'>"+genes[i]+"</a>,<br />";
      }
    }
    else {
      tip += ",<br /><strong>No genes from query taxon.</strong>,<br />";
    }
  }

  tip += ",<br />Click <a href='http://phylofacts.berkeley.edu/phog/netscope/"+node;

  if (taxon != "None") {
    tip += "_"+taxon;
  }

  tip += "/'>here</a> to re-center on this node...";
  tip_title = ""+node;
  Tip(tip, STICKY, true, CLOSEBTN, true, FONTSIZE, '9pt', OFFSETX, 4, OFFSETY, 4, TITLE, tip_title);

}

function edgeClick(edge, node1, node2) {
  if ("" != edge && "x" != edge) {
    window.open(edge);
  }
}

function edgeTip(edge, node1, node2) {
  if ("" != edges[edge][2]) {
    str = "Observed interaction between proteins "+node1+" and "+node2+".";
    if ("x" != edges[edge][2]) {
      str += "<ul><li>Click edge to open DIP source page...</li></ul>";
    }
  }
  else {
    str = "Predicted interaction between proteins "+node1+" and "+node2+".<br />";
    if (edges[edge][3][0][0] == edges[edge][3][0][1]) {
      str += "Based on in-species paralog relationship."
    }
    else {
      str += "Based on observed interactions between:<ul>";
      for (s in edges[edge][3]) {
	str += "<li>"+edges[edge][3][s][0]+" and "+edges[edge][3][s][1]+"</li>";
      }
    }
    str += "</ul>"
  }
  bpgToolTip(str);
}

function drawGraph(edge_list, node_coords, filename, taxon_scientific, taxon_common, query_id, ortholog_type) {

  taxon_scientific_name = taxon_scientific;
  taxon_common_name = taxon_common;
  ortho_type = ortholog_type;

  file_suffix = filename.substring(13);

  nodes = node_coords['nodes'];
  var node_count = nodes.length;
  var edge_count = edge_list.length;
  edges = edge_list;

  //document.write("NODE COUNT: "+node_count+" EDGE COUNT: "+edge_list.length+"<br>")
  
  //alert(query);
  //for (var i = 0; i < edge_count; i++) {
  //  for (j in edge_list[i]) {
  //    document.write(edge_list[i][j]);
  //    document.write(" | ");
  //  }
  //  document.write("<br>---<br>");
  //}
  
  
  //for (var i = 0; i < node_count; i++) {
  //  document.write(nodes[i]); 
  //  document.write(",<br />");
  //}
  
  document.write("<center>");
  //document.write("<strong>QUERY SEQUENCE ID:</strong> "+query_id+",  <strong>TAXON:</strong> "+taxon_scientific+" ("+taxon_common+"),<br />,<br /><hr><br,<br />");
  document.write("<div style='position:relative' name='network_div'><img border=0 src=\"http://phylogenomics.berkeley.edu"+file_suffix+".dat.png\" usemap='#mapA' width='800' height='600' alt=\"\"/>");
  document.write("<map name='mapA' id='mapA'>");
  for (var i = 0; i < node_count; i++) {
    document.write("<area shape=\"circle\" coords=\""+parseInt(nodes[i][1])+","+parseInt(nodes[i][2])+",12\" href=\"javascript:void(0);\" language=\"JavaScript\" onmouseover=\"nodeTip('"+nodes[i][4]+"','"+nodes[i][3]+"');return true;\"  onclick=\"nodeClick('"+nodes[i][4]+","+nodes[i][3]+"');return true;\" title=\"\" onmouseout=\"UnTip();return true;\">");
  }
  
  nodes.sort();
  
  var x1 = 0;
  var x2 = 0;
  var x3 = 0;
  var x4 = 0;
  var y1 = 0;
  var y2 = 0;
  var y3 = 0;
  var y4 = 0;
  
  cas = 0
  for (var i = 0; i < node_count; i++) {
    for (var j = 0; j < i; j++) {
      // check for interaction and write if present.
      //document.write("try "+nodes[i][0]+" <--> "+nodes[j][0]+"<br>");
      for (var k = 0; k < edge_count; k++) {
	if ( ((nodes[i][0] == edge_list[k][0]) &&
	      (nodes[j][0] == edge_list[k][1])) ||
	     ((nodes[i][0] == edge_list[k][1]) &&
	      (nodes[j][0] == edge_list[k][0])) ) {

	  xi = parseInt(nodes[i][1]); // x coordinate of node i
	  yi = parseInt(nodes[i][2]); // y coordinate of node i
	  xj = parseInt(nodes[j][1]); // x coordinate of node j
	  yj = parseInt(nodes[j][2]); // y coordinate of node j

	  cas = 0
	  if (xi < xj) {
	    if (yi < yj) {
	      x1 = xi - 4;
	      y1 = yi + 4;
	      x2 = xi + 4;
	      y2 = yi - 4;
	      x3 = xj + 4;
	      y3 = yj - 4;
	      x4 = xj - 4;
	      y4 = yj + 4;
	      cas = 1
	    }
	    else {
	      x1 = xi - 4;
	      y1 = yi - 4;
	      x2 = xi + 4;
	      y2 = yi + 4;
	      x3 = xj + 4;
	      y3 = yj + 4;
	      x4 = xj - 4;
	      y4 = yj - 4;
	      cas = 2
	    }
	  }
	  else {
	    if (yi < yj) {
	      x1 = xi - 4;
	      y1 = yi - 4;
	      x2 = xi + 4;
	      y2 = yi + 4;
	      x3 = xj + 4;
	      y3 = yj + 4;
	      x4 = xj - 4;
	      y4 = yj - 4;
	      cas = 3
	    }
	    else {
	      x1 = xi - 4;
	      y1 = yi + 4;
	      x2 = xi + 4;
	      y2 = yi - 4;
	      x3 = xj + 4;
	      y3 = yj - 4;
	      x4 = xj - 4;
	      y4 = yj + 4;
	      cas = 4
	    }
	  }
	  
	  //document.write("coords: "+x1+","+x2+" "+y1+","+y2+"<br>");

	  document.write("<area shape='poly' coords='"+x1+","+y1+","+x2+","+y2+","+x3+","+y3+","+x4+","+y4+"' href=\"javascript:void(0);\" lang=\"JavaScript\" title=\"\" onmouseover=\"edgeTip("+k+",'"+nodes[i][0]+"','"+nodes[j][0]+"');return true;\" onclick=\"edgeClick('"+edge_list[k][2]+"','"+nodes[i]+"','"+nodes[j]+"');return true;\" onmouseout=\"UnTip(); return true;\" alt=\"\" />");
	  
	  //alert("case: "+cas+"\n"+nodes[i][0]+"    "+nodes[j][0]+"\n"+nodes[i][1]+", "+nodes[i][2]+"\n"+nodes[j][1]+", "+nodes[j][2]+"\n"+x1+", "+y1+"\n"+x2+", "+y2+"\n"+x3+", "+y3+"\n"+x4+","+y4);
	  
	  break;
	}
      }
    }
  }
  
  document.write("</map>");
  document.write("</div>");
  //document.write(node_count+"<br />");
  //document.write(nodes+"<br />");
  //document.write("<br /><br />");
  document.write("</center>");
  
  //for (var i = 0; i < node_count; i++) {
  //document.write(nodes[i]+"<br />");
  //}
  
}

function drawPHOGGraph(node_coords, filename, taxon_id) {

  if (1 > edge_count) {
     return;
  }

  if (taxon_id) {
    taxon = taxon_id;
  }
  else {
    taxon = "None";
  }

  file_suffix = filename.substring(13);

  nodes = node_coords['nodes'];
  var node_count = nodes.length;
  
  /*
  //alert(query);
  for (var i = 0; i < edge_count; i++) {
    document.write(edge_list[i][0][0]);
    document.write(",");
    document.write(node_list[edge_list[i][0][0]]['description']);
    document.write(",");
    document.write(edge_list[i][0][1]);
    document.write(",");
    document.write(node_list[edge_list[i][0][1]]['description']);
    document.write(",0.95");
    document.write(",<br />");
  }
  */
  /*
  for (var i = 0; i < node_count; i++) {
    var genes = new Array();
    genes = node_list[nodes[i][0]]['genes_from_taxon'];
    document.write(node_list[nodes[i][0]]); 
    document.write(",<br />");
    if (genes) {
      document.write("number of taxon specific genes in "+nodes[i][0]+" = "+genes.length+",<br />");
      for (var j = 0; j < genes.length; j++) {
	document.write(genes[j]);
	document.write(",<br />");
      }
    }
  }
  */
  document.write("<center>");
  document.write("<div style='position:relative' name='network_div'><img src=\"http://phylogenomics.berkeley.edu"+file_suffix+".dat.png\" usemap='#mapA' width='696' height='648' alt=\"\"/>");
  document.write("<map name='mapA' id='mapA'>");
  for (var i = 0; i < node_count; i++) {
    document.write("<area shape=\"circle\" coords=\""+parseInt(nodes[i][1])+","+parseInt(nodes[i][2])+",12\" href=\"javascript:void(0);\" language=\"JavaScript\" onmouseover=\"nodeTip('"+nodes[i][0]+"');return true;\"  onclick=\"nodeClick('"+nodes[i][0]+"');return true;\" title=\"\" onmouseout=\"UnTip();return true;\">");
  }
  
  nodes.sort();
  
  var x1 = 0;
  var x2 = 0;
  var x3 = 0;
  var x4 = 0;
  var y1 = 0;
  var y2 = 0;
  var y3 = 0;
  var y4 = 0;
  
  for (var i = 0; i < node_count; i++) {
    // Fill in the entries where interac
    for (var j = 0; j < i; j++) {
      // check for interaction and write if present.
      for (var k = 0; k < edge_count; k++) {
	if ( ((nodes[i][0] == edge_list[k][0][0]) &&
	      (nodes[j][0] == edge_list[k][0][1])) ||
	     ((nodes[i][0] == edge_list[k][0][1]) &&
	      (nodes[j][0] == edge_list[k][0][0])) ) {
	  
	  if (nodes[i][1] < nodes[j][1]) {
	    if (nodes[i][2] < nodes[j][2]) {
	      x1 = parseInt(nodes[i][1]) + 4;
	      y1 = parseInt(nodes[i][2]) + 4;
	      x2 = parseInt(nodes[i][1]) - 4;
	      y2 = parseInt(nodes[i][2]) - 4;
	      x3 = parseInt(nodes[j][1]) - 4;
	      y3 = parseInt(nodes[j][2]) - 4;
	      x4 = parseInt(nodes[j][1]) + 4;
	      y4 = parseInt(nodes[j][2]) + 4;
	    }
	    else {
	      x1 = parseInt(nodes[i][1]) - 4;
	      y1 = parseInt(nodes[i][2]) + 4;
	      x2 = parseInt(nodes[i][1]) + 4;
	      y2 = parseInt(nodes[i][2]) - 4;
	      x3 = parseInt(nodes[j][1]) + 4;
	      y3 = parseInt(nodes[j][2]) - 4;
	      x4 = parseInt(nodes[j][1]) - 4;
	      y4 = parseInt(nodes[j][2]) + 4;
	    }
	  }
	  else {
	    if (nodes[i][2] < nodes[j][2]) {
	      x1 = parseInt(nodes[i][1]) - 4;
	      y1 = parseInt(nodes[i][2]) + 4;
	      x2 = parseInt(nodes[i][1]) + 4;
	      y2 = parseInt(nodes[i][2]) - 4;
	      x3 = parseInt(nodes[j][1]) + 4;
	      y3 = parseInt(nodes[j][2]) - 4;
	      x4 = parseInt(nodes[j][1]) - 4;
	      y4 = parseInt(nodes[j][2]) + 4;
	    }
	    else {
	      x1 = parseInt(nodes[i][1]) + 4;
	      y1 = parseInt(nodes[i][2]) + 4;
	      x2 = parseInt(nodes[i][1]) - 4;
	      y2 = parseInt(nodes[i][2]) - 4;
	      x3 = parseInt(nodes[j][1]) - 4;
	      y3 = parseInt(nodes[j][2]) - 4;
	      x4 = parseInt(nodes[j][1]) + 4;
	      y4 = parseInt(nodes[j][2]) + 4;
	    }
	  }
	  
	  document.write("<area shape='poly' coords='"+x1+","+y1+","+x2+","+y2+","+x3+","+y3+","+x4+","+y4+"' href=\"javascript:void(0);\" lang=\"JavaScript\" title=\"\" onmouseover=\"edgeTip('"+nodes[i][0]+"','"+nodes[j][0]+"');return true;\" onclick=\"edgeClick('"+nodes[i][0]+"','"+nodes[j][0]+"');return true;\" onmouseout=\"UnTip(); return true;\" alt=\"\" />,<br />");
	  
	  break;
	}
      }
    }
  }
  

  document.write("</map>");
  document.write("</div>");
  
  //document.write(node_count+"<br />");
  //document.write(nodes+"<br />");
  document.write("</center>");
  
  //for (var i = 0; i < node_count; i++) {
  //document.write(nodes[i]+"<br />");
  //}
  
}



function htmlEntityDecode(str) {

  //document.write(str);

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


/************* 
 * Input:  A list of edges representing pairs of interacting PHOGs.           
 *         Each edge is represented by it's pair of enpoints (PHOGs)          
 *         and the type of interaction.                                       
 *                                                                            
 * Output: A simple html table with the unique PHOGs in the list labelling    
 *         both the rows and columns.  For each edge, the type of interaction 
 *         is identified in the corresponding table entry.                    
 *                                                                            
 * Assumptions: A single interactio type.                                     
 *************/

function drawJsonHyperObj(aJsonHyperObj) {

  edge_list = aJsonHyperObj['edges'];
  edge_count = edge_list.length;

  if (1 > edge_count) {
     document.write("<h1>Either no interactions were found or too many were found.</h1>");
     return;
  }

  node_list = aJsonHyperObj['nodes'];
  node_count = node_list.length;
  query = aJsonHyperObj['query'];
  //document.write(query);
  //document.write("<br /><br />");
  //document.write(node_list);
  //document.write("<br /><br />");
  //document.write(edge_list);
  //document.write("<br /><br />");

  //document.write("EDGE COUNT: ");
  //document.write(edge_count);
  //document.write("<br /><br />");

  // Put all endoints (PHOGs) into a list, then reduce it to unique PHOGs
  var redundant_phog_count = edge_count * 2;
  var redundant_phog_list = new Array(redundant_phog_count);
  for (var i = 0; i < edge_count; i++) {
    redundant_phog_list[i] = edge_list[i][0][0];
  }
  for (var i = 0; i < edge_count; i++) {
    redundant_phog_list[i + edge_count] = edge_list[i][0][1];
  }

  // Write out the list of input edges (for debugging purposes).
  // Sorted, for easier viewing.
  edge_list.sort();
  //document.write("EDGES:<br />");
  //for (var i = 0; i < edge_count; i++) {
  //document.write(i + " ");
  //document.write(edge_list[i]);
  //document.write("<br />");
  //}
  //document.write("<br />");

  // Sort the redundant list for easy identification of duplicates.
  redundant_phog_list.sort();

  // Scan the redundant list and count unique PHOGs.
  // This is for allocation purposes and may not be needed (I'm new to js!)
  //var unique_phog_count = 1;
  for (var i = 1; i < redundant_phog_count; i++) {
    //document.write(i + "  ");
    //document.write(redundant_phog_list[i]);
    //document.write("<br />");
    if (redundant_phog_list[i] != redundant_phog_list[i-1]) {
      unique_phog_count++;
    }
  }
  //document.write("<br />");

  //document.write("NUMBER OF UNIQUE PHOGS: ");
  //document.write(unique_phog_count);
  //document.write("<br />");

  // Add one entry in the PHOG list for each unique PHOG.
  var unique_phog_count = 1;
  unique_phog_list[0] = redundant_phog_list[0];
  for (var i = 1; i < redundant_phog_count; i++) {
    if (redundant_phog_list[i] != redundant_phog_list[i-1]) {
      //unique_phog_list[unique_phog_count] = redundant_phog_list[i];
      unique_phog_list.push(redundant_phog_list[i]);
      unique_phog_count++;
    }
  }

  // Write debugging info.
  //document.write("NUMBER OF UNIQUE PHOGS: ");
  //document.write(unique_phog_count);
  //document.write("<br /><br />UNIQUE PHOGS:<br />");
  //for (var i = 0; i < unique_phog_count; i++) {
  //document.write(unique_phog_list[i]);
  //document.write("<br />");
  //}
  //document.write("<br />");

  document.write("<br />");


/*
  // Write the interaction table.
  document.write("<p>");
  document.write("<br />");
  document.write("<center>");
  document.write("<table border=\"1\" cellpadding=\"5\">");
  document.write("<tr>");
  document.write("<td>");
  document.write("</td>");

  for (var i = 0; i < unique_phog_count; i++) {
    document.write("<td align=\"center\">");
    document.write(i);
    document.write("</td>");
  }
  document.write("</tr>");

  // A column header for each PHOG.
  document.write("<tr><td></td>");
  for (var i = 0; i < unique_phog_count; i++) {
    document.write("<td align=\"center\">");
    var namelen = unique_phog_list[i].length;
    for (var j = 0; j < namelen; j++) {
      document.write(unique_phog_list[i][j]);
      document.write(",<br />");
    }
    document.write("</td>");
  }
  document.write("</tr>");


  // A row for each PHOG.
  for (var i = 0; i < unique_phog_count; i++){
    document.write("<tr>");
    document.write("<td align=\"right\">");
    document.write(unique_phog_list[i]);
    document.write("</td>");


    // Fill in the entries where interactions are present.
    for (var j = 0; j < i; j++) {
      document.write("<td align=\"center\">");

      // check for interaction and write type to table entry if present.
      for (var k = 0; k < edge_count; k++) {
	if ( ((unique_phog_list[i] == edge_list[k][0][0]) && (unique_phog_list[j] == edge_list[k][0][1])) || 
	     ((unique_phog_list[i] == edge_list[k][0][1]) && (unique_phog_list[j] == edge_list[k][0][0])) ) {
	  document.write(edge_list[k][1]);
	  break;
	}
      }

      document.write("</td>");
    }
    // Blank out the upper-right portion of the table for easier reading (I think).
    for (var j = i; j < unique_phog_count; j++) {
      document.write("<td bgcolor=\"grey\"></td>");
    }
    document.write("</tr>");
  }

  document.write("</table>");

  //document.write("<br />");
  document.write("<br />");
  document.write("<br />");
  document.write("</center>");
  document.write("</p>");
*/

  //graphinit(edge_list, edge_count, unique_phog_list, unique_phog_count);

  //appletinit(edge_list, edge_count, unique_phog_list, unique_phog_count);

  // return?
}
