function Label(constructorData) {
    /*
        constructorData should contain:

        x and y coordinates - these represent the the left side midpoint of a rectangle containing the label, 
            if these are not passed, the label is rendered at x = half the width of the paper and y= 20 pixels.

        If it is a leaf label, it should contain an object called leafLabel.
        otherwise it should contain an object called collapsedLabel.

        Inside leafLabel, 
        

    */
    if (constructorData) {
        if (constructorData.coords) {
            this.x = constructorData.coords.x;
            this.y = constructorData.coords.y;
        }
        else {this.x = Math.floor(phylogram.paperWidth/2); this.y = 20;}

        this.label = phylogram.paper.set();
        
        if (constructorData.taxonID) {
            this.taxonID = parseInt(constructorData.taxonID);
        }
        if (constructorData.numLeaves) {
            this.isCollapsedLabel = true;
            this.collapsedDescriptionObject = phylogram.paper.text(this.x, (this.y - phylogram.options.leafLineHeight),
                constructorData.numLeaves.toString() + " sequences from ").attr({'text-anchor': 'start'});
            this.label.push(this.collapsedDescriptionObject);
            this.taxonScientificObject = constructorData.taxonScientific ?
                phylogram.paper.text(this.label.getBBox().x2+1, (this.y - phylogram.options.leafLineHeight),
                    constructorData.taxonScientific) : phylogram.paper.text(this.label.getBBox().x2+1,
                    this.y - phylogram.options.leafLineHeight, "Unknown Taxa");
            this.taxonScientificObject.attr({'font-style':'italic', 'cursor':'default', 'text-anchor':'start'});
            if ((this.taxonID) && (this.taxonID > 2)) { //should this be something else?
                this.taxonScientificObject.attr({'href':'http://www.uniprot.org/taxonomy/' + this.taxonID.toString(),
                                                'cursor':'pointer',
                    'stroke':phylogram.options.linkColor});
            }
            this.label.push(this.taxonScientificObject);
        }
        else {
            this.isCollapsedLabel = false;
            this.taxonScientificObject = constructorData.taxonScientific ?
                phylogram.paper.text(this.x, (this.y - phylogram.options.leafLineHeight),
                    constructorData.taxonScientific) : phylogram.paper.text(this.x,
                    this.y - phylogram.options.leafLineHeight, "Unknown Taxon");
            this.taxonScientificObject.attr({'font-style':'italic', 'cursor':'default', 'text-anchor':'start'});
            if ((this.taxonID) && (this.taxonID > 2)) { //should this be something else?
                this.taxonScientificObject.attr({'href':'http://www.uniprot.org/taxonomy/' + this.taxonID.toString(),
                                                'cursor':'pointer',
                    'stroke':phylogram.options.linkColor});
            }
            this.label.push(this.taxonScientificObject);
        }
        this.taxonCommonObject = constructorData.taxonCommon ?
            phylogram.paper.text(this.label.getBBox().x2+1, 
                                    (this.y - phylogram.options.leafLineHeight), 
                                    " (" + constructorData.taxonCommon + ")").attr({'cursor':'default', 'text-anchor':'start'}) : null; 
        if (this.taxonCommonObject) {
            if ((this.taxonID) && (this.taxonID > 2)) {
                this.taxonCommonObject.attr({'href':'http://www.uniprot.org/taxonomy/' + this.taxonID.toString(),
                            'cursor':'pointer','stroke':phylogram.options.linkColor});
            }
            this.label.push(this.taxonCommonObject);
        }
        if (!(this.isCollapsedLabel)) {
            this.uniprotIDObject = constructorData.uniprotID ? 
            phylogram.paper.text(this.label.getBBox().x2+3, this.y-phylogram.options.leafLineHeight, 
                constructorData.uniprotID).attr({'text-anchor':'start', 'font-style':'bold','cursor':'pointer','href':'/phylofacts/sequence/UniProt/' + constructorData.uniprotID,'stroke':phylogram.options.linkColor}) : null;
            if (this.uniprotIDObject) {
                this.label.push(this.uniprotIDObject);
            }
        }
        this.swissProtIcon = constructorData.swissProt ?
            phylogram.paper.image('/static/img/icons/icon_swiss_flag_12.png', this.label.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                12, 12).attr({'text-anchor':'start'}) : null;
        if (this.swissProtIcon) {this.label.push(this.swissProtIcon);}
        this.GOIcon = constructorData.go ?
            phylogram.paper.image('/static/img/icons/icon_flask_12.png', this.label.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                9, 12).attr({'text-anchor':'start'}) : null;
        if (this.GOIcon) {this.label.push(this.GOIcon);}
        this.structureIcon = constructorData.structure ?
            phylogram.paper.image('/static/img/icons/icon_pdb_12.png', this.label.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                12, 12).attr({'text-anchor':'start'}) : null;
        if (this.structureIcon) {this.label.push(this.structureIcon);}
        this.literatureIcon = constructorData.literature ?
            phylogram.paper.image('/static/img/icons/icon_literature_12.png', this.label.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                9, 12).attr({'text-anchor':'start'}) : null;
        if (this.literatureIcon) {this.label.push(this.literatureIcon);}
        this.uniprotDescriptionObject = phylogram.paper.text(this.x, (this.y + phylogram.options.leafLineHeight), constructorData.uniprotDescription).attr({'text-anchor':'start'});
        if (this.isCollapsedLabel) {this.uniprotDescriptionObject.attr({'stroke':'blue'});}
        this.label.push(this.uniprotDescriptionObject);            
    }

    this.redrawIcons = function(icons) {
        // remove the current icons
        if (this.swissProtIcon) {
            this.label.exclude(this.swissProtIcon);
            this.swissProtIcon.remove();
            this.swissProtIcon = null;
        }
        if (this.GOIcon) {
            this.label.exclude(this.GOIcon);
            this.GOIcon.remove();
            this.GOIcon = null;
        }
        if (this.structureIcon) {
            this.label.exclude(this.structureIcon);
            this.structureIcon.remove();
            this.structureIcon = null;
        }
        if (this.literatureIcon) {
            this.label.exclude(this.literatureIcon);
            this.literatureIcon.remove();
            this.literatureIcon = null;
        }
        // draw the new icons
        var s = phylogram.paper.set();

        if (this.uniprotIDObject) {s.push(this.uniprotIDObject);}
        else if (this.taxonCommonObject) {s.push(this.taxonCommonObject);}
        else {s.push(this.taxonScientificObject);}

        if (icons.swissProt) {
            this.swissProtIcon = phylogram.paper.image('/static/img/icons/icon_swiss_flag_12.png', s.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                12, 12).attr({});
            s.push(this.swissProtIcon);
            this.label.push(this.swissProtIcon);
        }
        if (icons.go) {
            this.GOIcon = phylogram.paper.image('/static/img/icons/icon_flask_12.png', s.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                9, 12).attr({});
            s.push(this.GOIcon);
            this.label.push(this.GOIcon);
        }
        if (icons.structure) {
            this.structureIcon = phylogram.paper.image('/static/img/icons/icon_pdb_12.png', s.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                12, 12).attr({});
            s.push(this.structureIcon);
            this.label.push(this.structureIcon);
        }
        if (icons.literature) {
            this.literatureIcon = phylogram.paper.image('/static/img/icons/icon_literature_12.png', s.getBBox().x2+2, this.y - phylogram.options.leafLineHeight - 7,
                9, 12).attr({});
            this.label.push(this.literatureIcon);
        }
        s.clear();
    }

    this.remove = function() {
        if (this.collapsedDescriptionObject) {this.collapsedDescriptionObject.remove(); this.collapsedDescriptionObject = null;}
        if (this.uniprotIDObject) {this.uniprotIDObject.remove(); this.uniprotIDObject = null;}
        if (this.uniprotDescriptionObject) {this.uniprotDescriptionObject.remove(); this.uniprotDescriptionObject = null;}
        if (this.taxonScientificObject) {this.taxonScientificObject.remove(); this.taxonScientificObject = null;}
        if (this.taxonCommonObject) {this.taxonCommonObject.remove(); this.taxonCommonObject = null;}
        if (this.swissProtIcon) {this.swissProtIcon.remove(); this.swissProtIcon = null;}
        if (this.GOIcon) {this.GOIcon.remove(); this.GOIcon = null;}
        if (this.literatureIcon) {this.literatureIcon.remove(); this.literatureIcon = null;}
        if (this.structureIcon) {this.structureIcon.remove(); this.structureIcon = null;}
        if (this.label) {this.label.clear(); this.label = null;} 
    }
}
