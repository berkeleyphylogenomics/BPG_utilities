function Node(constructorData) {
    // defines the node object
    this.x = 0;
    this.y = 0;
    this.percentID = 0.0;
    this.isCollapsed = false;
    this.parentNode = null;
    this.children = [];
    this.branchLengthToRoot = 0;
    this.containedLeaves = 0;
    this.pathToParent = null;
    this.labelObject = null;
    this.nodeObject = null;
    this.glowObject = null;
    this.nodeObjectDefaultAttributes = constructorData.defaultAttr ? 
        constructorData.defaultAttr : {'cursor': 'pointer', 'fill-opacity':0, 'stroke-opacity':0,'r':6};
    this.isSelected = false;
    this.isHighlighted = false;
    var that=this;
    this.labelWidth = 1;   
    this.description = constructorData["name"];
    this.id = constructorData["nodeID"];
    this.branchLengthToParent = constructorData["branchLength"];
    this.hasGO = false;
    this.hasLiterature = false;
    this.hasStructure = false;
    this.hasSwissProt = false;
 
    if (constructorData["sequences"]) {
        this.uniprotAccession = constructorData["sequences"][0].accession ? 
            constructorData["sequences"][0].accession[0].Text : null;
        this.uniprotID = constructorData["sequences"][0].symbol ?
            constructorData["sequences"][0].symbol[0].Text : null;
        if (constructorData["sequences"][0].annotation) {
            for (var i=0; i<constructorData["sequences"][0].annotation.length; i++) {
                if (constructorData["sequences"][0].annotation[i].type == "go") {this.hasGO = true;}
                if (constructorData["sequences"][0].annotation[i].type == "literature") {this.hasLiterature = true;}
                if (constructorData["sequences"][0].annotation[i].type == "structure") {this.hasStructure = true;}
                if (constructorData["sequences"][0].annotation[i].type == "swissProt") {this.hasSwissProt = true;}
            }
        }
    }
    if (constructorData["taxonomy"]) {
        this.scientific_name = constructorData["taxonomy"][0].scientific_name[0].Text;
        this.common_name = constructorData["taxonomy"][0].common_name ? 
            constructorData["taxonomy"][0].common_name[0].Text : null;
        this.taxon_id = constructorData["taxonomy"][0].id ?
            parseInt(constructorData["taxonomy"][0].id[0].Text) : null;
    }
    if (constructorData["property"]) {
        // only property now is percent ID, set it.
        this.percentID = parseFloat(constructorData["property"][0].Text);
    }

    this.scaledX = function() {
        return Math.floor(this.x*phylogram.tree.horizontalScalingFactor+phylogram.options.leftPadding);
    }

    this.scaledY = function() {
        return Math.floor(this.y*phylogram.options.heightPerLeaf+phylogram.options.topPadding);
    }

    this.draw = function(attributes) {
        if (this.nodeObject) {this.nodeObject.remove();}
        if (attributes) {
            this.nodeObject = phylogram.paper.circle(this.scaledX(), this.scaledY(), 1).attr(attributes).toFront();
        }
        else {
            this.nodeObject = phylogram.paper.circle(this.scaledX(), this.scaledY(), 1).attr(this.nodeObjectDefaultAttributes).toFront();
        }
        if (this.isLeaf()) {
            
        }
        return this;
    }

    this.glow = function(attributes) {
        if (this.glowObject) {this.unglow();}
        var theseGlowAttr = attributes ? attributes : {};
        this.glowObject = this.nodeObject.glow(theseGlowAttr).toFront();
        return this;
    }

    this.unglow = function() {
        if (this.glowObject) {this.glowObject.remove(); this.glowObject = null;}
        return this;
    }

    this.setParent = function(parent) {this.parentNode = parent;}
    this.getParent = function() {return this.parentNode;}
    this.getChildren = function() {return this.children;}
    this.appendChild = function(child) {this.children.push(child);}
    this.isLeaf = function() {return (this.children.length === 0);}
    this.isRoot = function() {return !(this.parentNode);}
    this.hasSibling = function() {return (!(this.isRoot()) && (this.parentNode.children.length > 1));}
    this.makeThisLabel = function(coords) {
        var constructor = {};
        // guaranteed at this point to have a description
        constructor['uniprotDescription'] = this.description;
        if (this.uniprotID) {constructor['uniprotID'] = this.uniprotID;}
        if (this.scientific_name) {constructor['taxonScientific'] = this.scientific_name;}
        if (this.common_name) {constructor['taxonCommon'] = this.common_name;}
        if (this.taxon_id) {constructor['taxonID'] = this.taxon_id;}
        if (this.hasGO) {constructor['go'] = true;}
        if (this.hasLiterature) {constructor['literature'] = true;}
        if (this.hasStructure) {constructor['structure'] = true;}
        if (this.hasSwissProt) {constructor['swissProt'] = true;}
        if (!(this.isLeaf())) {constructor['numLeaves'] = this.containedLeaves;}
        if (coords) {constructor['coords'] = coords;}
        return constructor;
    }
}
