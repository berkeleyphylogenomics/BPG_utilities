function Tree(root) {
    // defines the tree object
    this.root = root;
    this.maximumBranchLength = 0;
    this.maximumDisplayedBranchLength = 0;
    var that = this;
    this.horizontalScalingFactor=0;
    this.leaves = [];
    this.displayed = [];
    this.longestBranchNode = null;
    this.longestBranchNodeLabelLength = 0;
    this.selectedNode = null;

    this.summarizeByDepth = function(depth) {
        function traverseTree(node,d) {
            if (d == depth) {
                node.isCollapsed = true;
            }
            else {
                node.isCollapsed = false;
                _.each(node.children, function(child) {
                    traverseTree(child, d+1);
                });
            }
        }
        traverseTree(this.root, 0);
    }

    this.expandAll = function() {
        function traverseTree(node) {
            node.isCollapsed = false;
            _.each(node.children, function(child) {
                traverseTree(child);
            });
        }
        traverseTree(this.root);
    }

    this.colorClade = function(node, attr) {
        function traverseTree(node) {
            _.each(node.children, function(child) {
                child.nodeObjectDefaultAttributes = attr.nodeAttr;
                if (!(child.isHighlighted)) {child.nodeObject.attr(attr.nodeAttr).toFront();}
                child.pathToParent.pathObjectDefaultAttributes = attr.pathAttr;
                if (!(child.pathToParent.isHighlighted)) {child.pathToParent.pathObject.attr(attr.pathAttr);}
                traverseTree(child);  
            });
        }  
        node.nodeObjectDefaultAttributes = attr.nodeAttr;
        if (!(node.isHighlighted)) {node.nodeObject.attr(attr.nodeAttr);}
        traverseTree(node);
    }

    this.colorByPercentID = function(threshold) {
        function traverseTree(node) {
            _.each(node.children, function(child) {
                if (child.percentID < threshold) {
                    child.nodeObjectDefaultAttributes = {'cursor': 'pointer', 'fill-opacity':0, 'stroke-opacity':0,'r':6};
                    if (!(child.isHighlighted)) {child.nodeObject.attr(child.nodeObjectDefaultAttributes).toFront();}
                    child.pathToParent.pathObjectDefaultAttributes = {'stroke-width':2, 'stroke':'black'};
                    if (!(child.pathToParent.isHighlighted)) {child.pathToParent.pathObject.attr(child.pathToParent.pathObjectDefaultAttributes);}
                    traverseTree(child);
                }
                else {
                    that.colorClade(child, {nodeAttr: {'cursor':'pointer','fill-opacity':0,'stroke-opacity':0,'r':6}, pathAttr: {'stroke-width':2,'stroke':phylogram.options.colorSet[colorIndex % phylogram.options.colorSet.length]}});
                    colorIndex++;
                }
            });
        }
        var colorIndex = 0;
        if (this.root.percentID > threshold) {
            if (threshold != 0) {this.colorClade(this.root, {nodeAttr: {'cursor':'pointer','fill-opacity':0,'stroke-opacity':0,'r':6}, pathAttr: {'stroke-width':2,'stroke':phylogram.options.colorSet[colorIndex % phylogram.options.colorSet.length]}});}
            else {this.colorClade(this.root, {nodeAttr: {'cursor':'pointer','fill-opacity':0,'stroke-opacity':0,'r':6}, pathAttr: {'stroke-width':2,'stroke':'black'}});}
        }
        else {traverseTree(this.root);}        
    }


    this.highlightPathToRoot = function(node, attr) {
        function traverseTree(node) {
            if (!(node.isRoot())) {
                if (attr) {
                    if (!(node.isHighlighted)) {node.nodeObject.attr(attr.nodeAttr).toFront();}
                    if (!(node.pathToParent.isHighlighted)) {node.pathToParent.pathObject.attr(attr.pathAttr);}
                }
                else {
                    if (!(node.isHighlighted)) {node.nodeObject.attr(node.nodeObjectDefaultAttributes).toFront();}
                    if (!(node.pathToParent.isHighlighted)) {node.pathToParent.pathObject.attr(node.pathToParent.pathObjectDefaultAttributes);}
                }
                traverseTree(node.parentNode);
            }
        }
        traverseTree(node);
    }

    this.glowPathToRoot = function(node, attr) {
        function traverseTree(node) {
            if (!(node.isRoot())) {
                if (attr) {
                    node.glow(attr.glowAttr);
                    node.pathToParent.glow(attr.glowAttr);
                    node.pathToParent.pathObject.attr(attr.pathAttr);
                    node.nodeObject.attr(attr.nodeAttr).toFront();
                    node.isHighlighted = true;
                    node.pathToParent.isHighlighted = true;
                }
                else {
                    node.unglow();
                    node.nodeObject.attr(node.nodeObjectDefaultAttributes);
                    node.pathToParent.unglow();
                    node.pathToParent.pathObject.attr(node.pathToParent.pathObjectDefaultAttributes);
                    node.isHighlighted = false;    
                    node.pathToParent.isHighlighted = false;
                }
                traverseTree(node.parentNode);
            }
        }
        traverseTree(node);
    }

    this.populateNodePositions = function() {
        function traverseTree(node) {
            if (!(node.isCollapsed || node.isLeaf())) {
                var dLeaves = 0;
                for(var i=0;i<node.children.length;i++) {
                    dLeaves+=traverseTree(node.children[i]);
                }
                node.displayedLeaves = dLeaves;
                node.x = node.branchLengthToRoot;
                if (node.children.length>1) {
                    node.y = (node.children[0].y + (node.children[node.children.length-1].y - node.children[0].y)/2);
                }
                else {node.y = node.children[0].y;}
                return dLeaves;
            }
            else
            {
                if (node.branchLengthToRoot > mdbl) {mdbl=node.branchLengthToRoot;}
                node.x = node.branchLengthToRoot;
                node.y = nLeaves;
                nLeaves+=1;
                var thisScalingFactor = (phylogram.paperWidth - node.labelWidth - phylogram.options.leftPadding)/node.x;
                if (thisScalingFactor < minScalingFactor) {that.longestBranchNode = node; minScalingFactor = thisScalingFactor;}
                return 1;
            }
        }
        var nLeaves = 0;
        var mdbl = 0;
        var minScalingFactor = 99999999;  // This is bad.
        traverseTree(this.root);
        this.maximumDisplayedBranchLength = mdbl;
        this.horizontalScalingFactor = minScalingFactor;
    }

    //this.colorTree = function(type,
    this.draw = function(options) {
        function recursivelyDraw(node) {
            if (!(node.isCollapsed || node.isLeaf())) {
                _.each(node.children, function(child) {
                    var pathString = "M" + node.scaledX().toString() 
                            + "," + node.scaledY().toString() +
                            "V" + child.scaledY().toString() + 
                            "H" + child.scaledX().toString();
                    if ((child.pathToParent) && (child.pathToParent.pathObjectDefaultAttributes)) {
                        child.pathToParent = new Path({'pathString': pathString, defaultAttr: child.pathToParent.pathObjectDefaultAttributes});
                    }
                    else {
                        child.pathToParent = new Path({'pathString': pathString, defaultAttr: {'stroke':'black','stroke-width':2}});
                    }
                    child.draw();
                    recursivelyDraw(child);
                });
            }
            else {
                // Put the label in the right spot and draw the node.
                node.draw();
                node.labelObject = new Label(node.makeThisLabel({'x': node.nodeObject.getBBox().x2, 'y': node.scaledY()}));
            }
            if (node != that.root) {
                node.nodeObject.mouseover(function() {
                    that.highlightPathToRoot(node, {pathAttr: {'stroke':'red'},
                                            nodeAttr:{'stroke':'red', 'fill':'red', 
                                            'stroke-opacity':1,'fill-opacity':1}});
                });
                node.nodeObject.mouseout(function() {
                    that.highlightPathToRoot(node);
                });
                node.nodeObject.click(function() {
                    if (node.isSelected) {
                        that.glowPathToRoot(node);
                        that.selectedNode = null;
                        if (node.isCollapsed) {node.isCollapsed = false;}
                        else {node.isCollapsed = true;}
                        node.isSelected = false;
                        that.draw();
                    }
                    else {
                        if (that.selectedNode) {
                            that.glowPathToRoot(that.selectedNode);
                            that.selectedNode.isSelected = false;
                        }
                        that.glowPathToRoot(node, {glowAttr: {'color':'red','width':6}, 
                            nodeAttr:{'stroke':'red','fill':'red','stroke-opacity':1,'fill-opacity':1,'r':6},
                            pathAttr: {'stroke':'red'}});
                        node.isSelected = true;
                        that.selectedNode = node;
                    }
                });
            }
        }
        if (options) {
            var refresh = options.refresh ? options.refresh : false;
        }
        // 15 to eliminate horizontal scrollbar
        phylogram.paperWidth = $j("#"+phylogram.options.div).width()-phylogram.options.scrollBarWidth;
        if (!(refresh)) {this.populateNodePositions();}
        else {this.horizontalScalingFactor = (phylogram.paperWidth - this.longestBranchNode.labelWidth - phylogram.options.leftPadding)/this.longestBranchNode.branchLengthToRoot;}
        phylogram.paperHeight = this.root.displayedLeaves*phylogram.options.heightPerLeaf+phylogram.options.topPadding;
        phylogram.paper.clear();
        phylogram.paper.setSize(phylogram.paperWidth, phylogram.paperHeight);
        recursivelyDraw(this.root);
        this.root.draw();
    }

    this.initializeTree = function() {
        // set paper to Raphael object, we will get the actualy canvas size later and set it.
        return this;
    }
}
