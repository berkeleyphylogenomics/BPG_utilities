function Path(constructorData) {
    this.pathObjectDefaultAttributes = constructorData.defaultAttr ? constructorData.defaultAttr : {};
    this.pathObject = constructorData.pathString ? phylogram.paper.path(constructorData.pathString).attr(this.pathObjectDefaultAttributes) : null;
    this.isHighlighted = false;

    this.setDefaultAttributes = function(pathAttr) {
        this.pathObjectDefaultAttributes = pathAttr ? pathAttr : {};
        return this; 
    }

    this.glow = function(attributes) {
        if (this.glowObject) {this.unglow();}
        var theseGlowAttr = attributes ? attributes : {};
        this.glowObject = this.pathObject.glow(theseGlowAttr);
        return this;
    }

    this.unglow = function() {
        if (this.glowObject) {this.glowObject.remove(); this.glowObject = null;}
        return this;
    }
}
