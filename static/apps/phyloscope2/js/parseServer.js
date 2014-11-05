function parseJSONtree(root) {
    // parses a json representation of the phyloxml
    function parseCladeData(clade) {
        var branchLength = Math.max(clade["branch_length"] || 1e-10, 1e-10);
        var name = clade.name ? clade.name[0].Text : "Unknown description";
        var cladeData = {
            "branchLength" : branchLength,
            "name" : name,
            "nodeID" : parseInt(clade.node_id[0].Text)
        };
        if (clade.sequence) {
            cladeData["sequences"] = clade.sequence;
        }
        if (clade.property) {
            cladeData["property"] = clade.property;
        }
        if (clade.taxonomy) {
            cladeData["taxonomy"] = clade.taxonomy;
        }
        return cladeData;
    }

    function parseClade(clade, parentNode) {
        var node = new Node(parseCladeData(clade));
        var cLeaves = 0;
        node.parentNode = parentNode;
        node.branchLengthToRoot = node.branchLengthToParent + node.parentNode.branchLengthToRoot;
        if (node.branchLengthToRoot > maxBranchLength) {maxBranchLength = node.branchLengthToRoot;}
        _.each(clade.clade, function(child) { 
            var retObj = parseClade(child, node);
            node.appendChild(retObj.node);
            cLeaves += retObj.leaves;
        });
        node.labelObject = new Label(node.makeThisLabel({'x':300, 'y':100}));
        node.labelWidth = node.labelObject.label.getBBox().width
        node.minScalingFactor = (phylogram.paperWidth - node.labelWidth - phylogram.options.leftPadding)/node.x;
        if (node.children.length == 0) {
            tree.leaves.push(node);
            tree.displayed.push(node);
            return {'node': node, 'leaves': 1, 'hasSwissProt': node.hasSwissProt,
                'hasGO': node.hasGO, 'hasStructure': node.hasStructure, 'hasLiterature': node.hasLiterature};
        }
        else {
            var sp = false;
            var go = false;
            var s = false;
            var l = false;
            _.each(node.children, function(child) {
                if (child.hasSwissProt) {sp = true;}
                if (child.hasGO) {go = true;}
                if (child.hasStructure) {s = true;}
                if (child.hasLiterature) {l = true;}
            });
            node.hasSwissProt = sp;
            node.hasGO = go;
            node.hasStructure = s;
            node.hasLiterature = l;
            node.containedLeaves = cLeaves;
            return {'node': node, 'leaves': cLeaves, 'hasSwissProt': node.hasSwissProt,
                'hasGO': node.hasGO, 'hasStructure': node.hasStructure, 'hasLiterature': node.hasLiterature};
        }
    }

    var numLeaves = 0;
    var maxBranchLength = 0;
    var rootNode = new Node(parseCladeData(root));
    var tree = new Tree(rootNode);

    _.each(root.clade, function(child) {
        var retObj = parseClade(child, rootNode);
        rootNode.appendChild(retObj.node);
        numLeaves+=retObj.leaves;
    });
    if (rootNode.children.length == 0) {tree.leaves.push(rootNode);}
    tree.root.containedLeaves = numLeaves;
    tree.maximumBranchLength = maxBranchLength;
    return tree;
}
