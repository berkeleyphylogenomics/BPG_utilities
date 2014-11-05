Smits = {};Smits.Common = {
    nodeIdIncrement : 0,
    activeNode: 0,
    
    /* Rounds float to a defined number of decimal places */
    roundFloat : function(num, digits){
        var i = 0, 
            dec = 1;
        while(i < digits){
            dec *= 10;
            i++;
        }
        return Math.round(num*dec)/dec; 
    },
    
    /* Copies properties from one object to another */
    apply : function(obj, extObj){
        if (obj && typeof extObj == 'object') {
            for (var key in extObj) {
                obj[key] = extObj[key];
            }
        }
        return obj; 
    },
    
    addRaphEventHandler : function(el, eventType, fn, paramsObj){
        try{
            el[eventType](function(fn, paramsObj){
                return function(e,o){
                    var params = paramsObj;
                    params.e = e;
                    fn(params);
                };
            }(fn, paramsObj));
        } catch (err){} 
    },
    
    isInteger : function(s) {
        return !isNaN(parseInt(s));
    },

    isXMLSerializerAvailable : function(){
        if (typeof(XMLSerializer) == "function"){
            return true;
        } else {
            return false;
        }
    },
    
    createSvgEl : function (el, attr) {
        el = document.createElementNS("http://www.w3.org/2000/svg", el);            
        if (attr) {
            for (var key in attr) {
                if (attr.hasOwnProperty(key)) {
                    el.setAttribute(key, String(attr[key]));
                }
            }
        }   
        return el;  
    },
    
    createGradientEl : function(name, obj, coords){
        if(obj.type != "radialGradient") return false;
        
        var radialEl = Smits.Common.createSvgEl("radialGradient", {
            id: name, 
            gradientUnits:"userSpaceOnUse", 
            cx: coords[0], 
            cy: coords[1], 
            r: coords[2], 
            fx: coords[0], 
            fy: coords[1]
        });

        if(obj.stop){
            var stop = obj.stop;
            for(var i = 0; i < stop.length; i++){
                var stopObj = stop[i];
                if(stopObj['@attributes']){
                    radialEl.appendChild(Smits.Common.createSvgEl("stop", stopObj['@attributes']));
                } else {
                    if(stopObj['_attributes']) delete stopObj['_attributes'];
                    if(stopObj['_children']) delete stopObj['_children'];
                    if(stopObj['__proto__']) delete stopObj['__proto__'];
                    radialEl.appendChild(Smits.Common.createSvgEl("stop", stopObj));
                }
            }
        }
        
        return radialEl;
    },
    
    setCssStyle : function(selector, rule) {
        var stylesheet = document.styleSheets[0];
        if( stylesheet.addRule ){
            stylesheet.addRule(selector, rule);
        } else if( stylesheet.insertRule ){
            stylesheet.insertRule(selector + ' { ' + rule + ' }', stylesheet.cssRules.length);
        }
    }
};

/* AUTOMATIC APPLICATION */
(function(){
    // Mozilla-provided fix for browsers without indexOf
    if (!Array.prototype.indexOf)
    {
      Array.prototype.indexOf = function(elt /*, from*/)
      {
        var len = this.length >>> 0;

        var from = Number(arguments[1]) || 0;
        from = (from < 0)
             ? Math.ceil(from)
             : Math.floor(from);
        if (from < 0)
          from += len;

        for (; from < len; from++)
        {
          if (from in this &&
              this[from] === elt)
            return from;
        }
        return -1;
      };
    }   
})();Smits.PhyloCanvas = function(inputFormat, sDivId, canvasWidth, canvasHeight, type){
    var phylogram,
        divId,
        newickObject,
        svg,
        dataObject,
        availableTypes = ['phylogram', 'circular'];

    var refresh = function(){
        svg.clear();
        render( refresh );
    }
        
    var render = function( phyloCanvasObject ){
        var namedType = availableTypes[type];
        switch (namedType){
            case "circular":
                phylogram = new Smits.PhyloCanvas.Render.CircularPhylogram(
                    svg, 
                    dataObject,
                    phyloCanvasObject
                );              
            break;
            
            default:
                phylogram = new Smits.PhyloCanvas.Render.Phylogram(
                    svg,
                    dataObject,
                    phyloCanvasObject
                );              
        }
    };
    
        
    /* CONSTRUCTOR */
    // Process dataset -- assume newick format, else needs to provide format
    if(typeof inputFormat === "object"){
        if(inputFormat.xml){    // default xml format is phyloXML
            if(!inputFormat.fileSource){
                var xj = XMLObjectifier.textToXML(inputFormat.xml);             // assume we need to clean it up
            } else {
                var xj = inputFormat.xml;
            }
            xj = XMLObjectifier.xmlToJSON(xj);
            dataObject = new Smits.PhyloCanvas.PhyloxmlParse(xj);
        } else if(inputFormat.phyloxml){
            if(!inputFormat.fileSource){
                var xj = XMLObjectifier.textToXML(inputFormat.phyloxml);            // assume we need to clean it up
            } else {
                var xj = inputFormat.phyloxml;
            }
            xj = XMLObjectifier.xmlToJSON(xj);
            dataObject = new Smits.PhyloCanvas.PhyloxmlParse(xj);
        } else if(inputFormat.nexml){
            if(!inputFormat.fileSource){
                var xj = XMLObjectifier.textToXML(inputFormat.nexml);           // assume we need to clean it up
            } else {
                var xj = inputFormat.nexml;
            }
            xj = XMLObjectifier.xmlToJSON(xj);
            dataObject = new Smits.PhyloCanvas.NexmlParse(xj, inputFormat)
        } else if(inputFormat.json){
            dataObject = new Smits.PhyloCanvas.PhyloxmlParse(inputFormat.json);
        } else if(inputFormat.newick){
            dataObject = new Smits.PhyloCanvas.NewickParse(inputFormat.newick);
        } else if(inputFormat.nexmlJson){
            dataObject = new Smits.PhyloCanvas.NexmlJsonParse(inputFormat);             
        } else {
            alert('Please set the format of input data');
        }
    } else {
        dataObject = new Smits.PhyloCanvas.NewickParse(inputFormat);
    }

    // If we have a phylogram, we want to set the minimum height high enough for all of the leaves
    if ( type !== 'circular' ) {
        var numLeaves = dataObject.getRoot().getCountAllChildren();
        var labelHeight = Smits.PhyloCanvas.Render.Style.getStyle('bootstrap', 'text')["font-size"];
        var minHeight = 20 + numLeaves * labelHeight * 2;
        if (canvasHeight < minHeight) { 
            canvasHeight = minHeight;
        }
    }

    divId = sDivId;
    svg = new Smits.PhyloCanvas.Render.SVG( divId, canvasWidth, canvasHeight );
    
    /* FACTORY */
    if( availableTypes.indexOf(type) != -1 ){
        type = availableTypes.indexOf(type) || 0;
    } else {
        type = 0;
    }
    render( refresh );
        
    return {

        getDataObject : function(){
            return dataObject;
        },
        scale : function(multiplier){
            svg.svg.scale(multiplier);
        },
        getSvg : function(){
            return svg;
        },
        getPhylogram : function(){
            return phylogram;
        },
        getSvgSource : function(){
            if(Raphael.svg && Smits.Common.isXMLSerializerAvailable()){
                var serialize = new XMLSerializer();
                return serialize.serializeToString(svg.svg.canvas);
            } else {
                return false;
            }
        },
        switchType : function(){
            if( type == availableTypes.length - 1 ){
                type = 0;
            } else {
                type++;
            }
            
            svg.clear();
            render();
        },
        
        refresh :refresh,
        render: render
    }

};
Smits.PhyloCanvas.Node = function(o, parentInstance){
    /**
    * Node Class
    * Allows objects to be traversed across children
    *
    */
    
    // initiate object
    var id = Smits.Common.nodeIdIncrement += 1;
    var level = 0;
    var len = 0;
    var newickLen = 0;
    var name = '';
    var type = '';
    var chart = {};
    var img = [];
    var collapsed = false;
    var collapsedCap = false;
    var children = [];
    
    if(o) Smits.Common.apply(this, o);

    /* Cache Calculations */
    var _countAllChildren = false;
    var _countImmediateChildren = false;
    var _midBranchPosition = false;
    var _offspring = false;
    var _allChildrenCombinedMaxLength = false;
    
    if(parentInstance){
        parentInstance.children.push(this); 
    }   
    
    return {
        id: id,
        level: level,
        len: len,
        newickLen: newickLen,
        name: name, 
        type: type,
        chart: chart,
        img: img,
        parentInstance: parentInstance,
        collapsed: collapsed,
        children: children,
        collapsedCap: collapsedCap,
        
        /* 
         * Only returns children that are to be rendered (e.g. non-collapsed)
         */
        getRenderChildren : function(){
            if( _offspring !== false ) return _offspring;

            var offspring = [];
            for (var i = 0; i < children.length; i++) {
                if( !children[i].collapsed ){
                    offspring.push( children[i] );
                }
            }
            return offspring;
        },
        
        getAllChildrenCombinedMaxLength : function( isChild ){
            if( _allChildrenCombinedMaxLength !== false ) return _allChildrenCombinedMaxLength;
            var len = 0;
            var combinedLen = [0];
            for( var i = 0; i < children.length; i++ ){
                combinedLen.push( children[i].getAllChildrenCombinedMaxLength(true) );
            }
            
            len = Math.max.apply(Math, combinedLen);
            if( isChild ){
                len += this.len;
            }
            return len;
        },
        
        getImmediateChildrenCombinedMinLength : function(){
            var combinedLen = [];
            for( var i = 0; i < children.length; i++ ){
                combinedLen.push( children[i].getAllChildrenCombinedMaxLength(true) );
            }
            return Math.min.apply(Math, combinedLen);
        },
        
        getCountAllChildren : function(){
            if( _countAllChildren !== false ) return _countAllChildren;
            var nodeCount = 0;

            var offspring = this.getRenderChildren();
            for (var i = 0; i < offspring.length; i++) {
                if(Smits.Common.isInteger(i)){
                    var child = offspring[i];
                    if( child.getCountAllChildren && child.getCountAllChildren() > 0 ){
                        nodeCount += child.getCountAllChildren();
                    } else {
                        nodeCount ++;
                    }           
                }
            }
            _countAllChildren = nodeCount;
            return nodeCount;
        },
        
        getMidbranchPosition : function(firstBranch){
            if( _midBranchPosition !== false ) return _midBranchPosition;
            var y = [0,0];  // bounds
            
            var offspring = this.getRenderChildren();
            for (var i = 0; i < offspring.length; i++) {
                var child = offspring[i];
                if( child.getCountAllChildren() ){
                    if(i == 0 && firstBranch){
                        y[0] = child.getMidbranchPosition(true);                
                        y[1] += child.getCountAllChildren() - 1;    
                    } else if(i == 0){
                        y[0] = child.getMidbranchPosition();                
                        y[1] += child.getCountAllChildren();    
                    } else if (i == offspring.length - 1){
                        y[1] += child.getMidbranchPosition();
                    } else {
                        y[1] += child.getCountAllChildren();                
                    }
                } else {
                    if(i == 0 && firstBranch){
                        y[0] = 0;
                    } else if(i == 0){
                        y[0] = 1;
                        y[1] += 1;  
                    } else if (i == offspring.length - 1){
                        y[1] += 1;
                    } else {
                        y[1] += 1;
                    }
                }
            }
            
            _midBranchPosition = y[1] >= y[0] ? ((y[1] + y[0]) / 2) : y[0];
            return _midBranchPosition;
        },

        clearCache : function(){
            _countAllChildren = false;
            _countImmediateChildren = false;
            _midBranchPosition = false;
            _offspring = false;     
            
            for (var i = 0; i < children.length; i++) {
                children[i].clearCache();
            }           
        }
    
    
    }
};

Smits.PhyloCanvas.PhyloxmlParse = function(jsonString){

    var mLevel = 0,
    mNewickLen = 0,
    root,
    validate,
    barCharts = [], binaryCharts = [], integratedBinaryCharts = [],
        
    recursiveParse = function(clade, parentNode){
        var node = new Smits.PhyloCanvas.Node();
        if(parentNode){
            node.level = parentNode.level + 1;
            //node.parentNode = parentNode;
        }
        if(clade.clade && clade.clade.length){
            for(var i = 0; i < clade.clade.length; i++){
                var thisClade = clade.clade[i];
                node.children.push(recursiveParse(thisClade, node));
            }
        }
        if(clade.branch_length){    // Branches can be attributes or own element
            if(typeof clade.branch_length === 'object'){
                clade.branch_length = clade.branch_length[0].Text;
            }

            node.len = Smits.Common.roundFloat(clade.branch_length, 4);         // round to 4 decimal places
            if(node.len == 0){
                node.len = 0.0001;
            }           
        }
        if(clade.name){
            node.type = 'label';
            node.name = clade.name[0].Text;
            if(clade.name[0] && clade.name[0].style){
                node.style = clade.name[0].style;
            }
            if(clade.name[0] && clade.name[0].bgStyle){
                node.bgStyle = clade.name[0].bgStyle;
            }           
        } else if(clade.confidence){
            node.name = clade.confidence[0].Text;
        }

        /* Collect further info that might be used as a label */
        if (clade.sequence && clade.sequence[0] && clade.sequence[0].name && clade.sequence[0].name[0] && clade.sequence[0].name[0].Text){
            node.sequenceName = clade.sequence[0].name[0].Text;
        }
        if (clade.taxonomy && clade.taxonomy[0]){
            if(clade.taxonomy[0].scientific_name && clade.taxonomy[0].scientific_name[0] && clade.taxonomy[0].scientific_name[0].Text){
                node.taxonomyScientificName = clade.taxonomy[0].scientific_name[0].Text;
            }
            if (clade.taxonomy[0].common_name  && clade.taxonomy[0].common_name[0] && clade.taxonomy[0].common_name[0].Text){
                node.taxonomyCommonName = clade.taxonomy[0].common_name[0].Text;
            }
            if (clade.taxonomy[0].id  && clade.taxonomy[0].id[0] && clade.taxonomy[0].id[0].Text){
                node.taxonomyID = clade.taxonomy[0].id[0].Text;
            }
        }
        if (clade.sequence && clade.sequence[0] && clade.sequence[0].accession && clade.sequence[0].accession[0] && clade.sequence[0].accession[0].Text){
            node.sequenceAccession = clade.sequence[0].accession[0].Text;
        }
        if (clade.point ){
            node.LatLong = [clade.point[0].lat[0].Text, clade.point[0]['long'][0].Text];
        }       

        
        /* Prioritization of Label */
        if(!node.name){
            if(node.sequenceName){
                node.name = node.sequenceName;
            } else if (node.taxonomyScientificName){
                node.name = node.taxonomyScientificName;
            } else if (node.taxonomyCommonName){
                node.name = node.taxonomyCommonName;
            } else if (node.sequenceAccession){
                node.name = node.sequenceAccession;
            }
            if(node.name){  // if name is now set, type is 'label'
                node.type = 'label'; 
            }
        }
        
        if( clade.collapsed && clade.collapsed[0] && clade.collapsed[0].Text ){
            node.collapsed = true;
        }
        
        if(clade.annotation){
            node.annotations = clade.annotation;
            if(clade.annotation[0] && clade.annotation[0].desc && clade.annotation[0].desc[0] && clade.annotation[0].desc[0].Text){
                node.description = clade.annotation[0].desc[0].Text;
            }
//            if(clade.annotation[0] && clade.annotation[0].uri && clade.annotation[0].uri[0] && clade.annotation[0].uri[0].Text){
//                node.uri = clade.annotation[0].uri[0].Text;
//            }           
            if(clade.annotation[0] && clade.annotation[0].img){
                for(var i in clade.annotation[0].img){
                    if(Smits.Common.isInteger(i)){
                        node.img[i] = clade.annotation[0].img[i].Text;
                    }
                }
            }
        }
        if(clade.chart){
            if(clade.chart[0]){
                for(var i in clade.chart[0]){
                    if(i != 'Text' && i != '_children'){
                    node.chart[i] = clade.chart[0][i][0].Text;
                    }
                }
            }
        }

        if (clade.color){
            node.color = clade.color[0];
            console.log(node);
        }
        
        // Validation
        if(node && node.level > 1){
            if(!node.len){
                validate = 'Error. Please include Branch Lengths - we only draw rooted phylogenetic trees.';
            }
        }
            
        return node;
    },
    
    recursiveProcessRoot = function(node, parentNode){
        if(node.children && node.children.length){
            for( var i = 0; i < node.children.length; i++){
                var child = node.children[i];
                child.newickLen = Math.round( (child.len + node.newickLen) *10000)/10000;
                if(child.level > mLevel) mLevel = child.level;
                if(child.newickLen > mNewickLen) mNewickLen = child.newickLen;
                if(child.children.length > 0){
                    recursiveProcessRoot(child, node); 
                }               
            }
        }
        return node;
    },
    
    recursiveProcessParameters = function(parametersEl, treeType){
        for (var i in parametersEl){
            if(i != '_children' && i != 'Text'){
                if(i == 'rectangular' || i == 'circular'){
                    recursiveProcessParameters(parametersEl[i][0], i);
                } else {
                    if(!Smits.PhyloCanvas.Render.Parameters[i]) {  Smits.PhyloCanvas.Render.Parameters[i] = {}; };
                    Smits.PhyloCanvas.Render.Parameters.set(i, parametersEl[i][0].Text, treeType);
                }
            }
        }
        return;
    };

    /* CONSTRUCTOR */   
    if(jsonString.phylogeny && jsonString.phylogeny[0] && jsonString.phylogeny[0].clade){
        root = recursiveParse(jsonString.phylogeny[0].clade[0]);
    }
    
    if(jsonString.phylogeny && jsonString.phylogeny[0] && jsonString.phylogeny[0].render && jsonString.phylogeny[0].render[0]){
        var render = jsonString.phylogeny[0].render[0];
        
        // Custom Styles
        if(render && render.styles){
            var styles = render.styles[0];
            for (var i in styles){
                if(i != '_children' && i != 'Text'){
                    if(styles[i][0]['type'] && styles[i][0]['type'] == "radialGradient" && Raphael.svg){
                        // radialGradient only supported by SVG
                        styles[i][0]['name'] = i;
                        Smits.PhyloCanvas.Render.Style[i] = styles[i][0];
                        if(!Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList']) { Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'] = [] };
                        Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'].push(i); 
                    } else {
                        if(!Smits.PhyloCanvas.Render.Style[i]) {  Smits.PhyloCanvas.Render.Style[i] = {}; };
                        for(var j in styles[i][0]){
                            if(j != '_attributes' && j != '_children' && j != 'type'){
                                Smits.PhyloCanvas.Render.Style[i][j.replace('_', '-')] = styles[i][0][j];       // This is quite painful, as xml does not allow dashes
                            }
                        }
                    }
                    
                }
            }
        }
        
        // Custom Parameters
        if(render && render.parameters){
            recursiveProcessParameters(render.parameters[0]);
        }           
        
        // Charts
        if(render && render.charts){
            var charts = render.charts[0];
            for (var i in charts){
                if(i != '_children' && i != 'Text'){
                    for(var j in charts[i]){
                        if(charts[i][j].type == "binary"){
                            charts[i][j].chart = i;
                            binaryCharts.push(charts[i][j]);
                        } else if (charts[i][j].type == "integratedBinary"){
                            charts[i][j].chart = i;
                            integratedBinaryCharts.push(charts[i][j]);                              
                        } else if (charts[i][j].type == "bar"){
                            charts[i][j].chart = i;
                            barCharts.push(charts[i][j]);
                        }
                    }
                }
            }
        }           
        
    }
    mNewickLen = root.len;
    root = recursiveProcessRoot(root);
    
    return {
        barCharts : barCharts, 
        binaryCharts : binaryCharts, 
        integratedBinaryCharts : integratedBinaryCharts,
        
        getRoot : function(){
            return root;
        },
        getLevels : function(){
            return mLevel;
        },
        getNewickLen : function(){
            return mNewickLen;
        },      
        getValidate : function(){
            return validate;
        },
        getBinaryCharts : function(){
            return binaryCharts;
        },
        getIntegratedBinaryCharts : function(){
            return integratedBinaryCharts;
        },
        getBarCharts: function(){
            return barCharts;
        }
        
    }

};
Smits.PhyloCanvas.Render = {};Smits.PhyloCanvas.Render.Style = {
    /* Default Styles */
    
    line: {
        "stroke":       'rgb(0,0,0)',
        "stroke-width": 1
    },
    
    text: {
        "font-family":  'Verdana',
        "font-size":    12,
        "text-anchor":  'start'
    },
    
    path: {
        "stroke":       'rgb(0,0,0)',
        "stroke-width": 1   
    },
    
    connectedDash : {
        "stroke":           'rgb(200,200,200)',
        "stroke-dasharray": ". "
    },
    
    textSecantBg : {
        "fill":     '#EEE',
        "stroke":   '#DDD'
    },
    
    highlightedEdgeCircle : {
        "stroke":   'red'
    },
    
    barChart : {
        fill:       '#003300',
        stroke:     '#DDD'
    },
    
    getStyle : function(requestStyle, fallbackStyle){
        if(this[requestStyle]){
            return this[requestStyle];
        } else {
            return this[fallbackStyle];
        }
    
    }
};Smits.PhyloCanvas.Render.Parameters = {

    /* DEFAULT PARAMETERS */
    jsOverride: 0,              // If set, js will override chart's file setting
    
    /** Phylogram parameters are separated because they behave very differently **/
    
    /* Rectangular Phylogram */
    Rectangular : {
        bufferX         : 200,          // Reduces the available canvas space for tree branches, allowing
                                        // for more space for the textual/charting components
        paddingX        : 10,
        paddingY        : 20,
        bufferInnerLabels : 10,         // Pixels
        bufferOuterLabels : 5,          // Pixels
        minHeightBetweenLeaves : 10,    // Should probably set pretty low, as clipping may occur if it needs to be implemented      
        
        alignPadding    : 0,            // Pixels to push the labels out by - this extension should be 
                                        // compensated by an increase in bufferX too
        alignRight      : false,
        
        showScaleBar    : false         // (STRING,  e.g. "0.05") Shows a scale bar at the bottom of the tree
    },
    
    /* Circular Phylogram */
    Circular : {
        bufferRadius        : 0.33,     // Margins of Tree Circle
                                        // If > 1, it is in pixels
                                        // If < 1, it is a percentage of the full canvas size       
        bufferAngle         : 20,       // controls split size in circle        
        initStartAngle      : 160,      
        innerCircleRadius   : 0,
        minHeightBetweenLeaves : 5,

        /* Labels */
        bufferInnerLabels : 2,          // Pixels
        bufferOuterLabels : 5,          // Pixels
        showBootstrapValues: false
    },
    
    /* Charts */

        /* Binary Defaults */
        binaryChartBufferInner : 5, 
        binaryChartBufferSiblings : 0.01,
        binaryChartThickness : 15,
        binaryChartDisjointed : false,
            
        /* Bar Defaults */
        barChartBufferInner : 3,
        barChartHeight : 50,
        barChartWidth : 0.5,    // If > 1, it is in pixels
                                // If < 1, it is a percentage of the node width 
                        
        /* 
            Rollover Events 
                At minimum, the params object has the following properties:
                    .svg
                    .node
                    .x
                    .y
                    .textEl
        */
        mouseRollOver : function(params) {
            if(params.node.edgeCircleHighlight){
                params.node.edgeCircleHighlight.show();
            } else {
                var circleObject = params.svg.draw(
                    new Smits.PhyloCanvas.Render.Circle(
                        params.x, params.y, 5,
                        { attr: Smits.PhyloCanvas.Render.Style.highlightedEdgeCircle }
                    )
                );
                params.node.edgeCircleHighlight = circleObject[0];
            }                   
            params.textEl.attr(Smits.PhyloCanvas.Render.Style['text-highlight']);
        },
        mouseRollOut : function(params) {
            params.node.edgeCircleHighlight.hide();
            params.textEl.attr(Smits.PhyloCanvas.Render.Style['text']);
        },

    set : function(param, value, treeType){
        if(!this.jsOverride){
            if(treeType){
                if(treeType == 'circular'){             
                    this['Circular'][param] = parseFloat(value);
                } else if (treeType == 'rectangular'){
                    this['Rectangular'][param] = parseFloat(value);
                }
            } else {
                this[param] = parseFloat(value);
            }
        }
    }
};Smits.PhyloCanvas.Render.Line = function(){

    return function(x1, x2, y1, y2, params){
        /* Defaults */  
        this.type = 'line';
        this.attr = Smits.PhyloCanvas.Render.Style.line;
        
        this.x1 = x1;
        this.x2 = x2;
        this.y1 = y1;
        this.y2 = y2;
        
        if(params) {
            Smits.Common.apply(this, params);
            if(params.attr) this.attr = params.attr;
        }

    }
}();Smits.PhyloCanvas.Render.Text = function(){

    return function(x, y, text, params){
        /* Defaults */
        this.type = 'text';
        this.attr = Smits.PhyloCanvas.Render.Style.text;
        
        this.x = x;
        this.y = y;
        this.text = text;
        
        if(params) {
            Smits.Common.apply(this, params);
            if(params.attr) this.attr = params.attr;
        }
    }
}();Smits.PhyloCanvas.Render.Path = function(){
    var attr = Smits.PhyloCanvas.Render.Style.path;
    
    return function(path, params){
        /* Defaults */
        this.type = 'path';
        this.attr = Smits.PhyloCanvas.Render.Style.path;
        
        this.path = path;
        if(params) {
            Smits.Common.apply(this, params);
            if(params.attr) this.attr = params.attr;
        }

    }
}();Smits.PhyloCanvas.Render.Circle = function(){

    return function(x, y, radius, params){
        /* Defaults */  
        this.type = 'circle';
    
        this.x = x;
        this.y = y;
        this.radius = radius;
        
        if(params) {
            Smits.Common.apply(this, params);
            if(params.attr) this.attr = params.attr;
        }
    }
}();Smits.PhyloCanvas.Render.SVG = function(divId, canvasWidth, canvasHeight){
    var canvasSize,
        svg;
        

    /* CONSTRUCTOR */
    canvasSize = [canvasWidth, canvasHeight];
    
    svg = Raphael(divId, canvasSize[0], canvasSize[1]);

        
    return {
        divId: divId,
        svg: svg,
        canvasSize: canvasSize,
        
        clear : function(){
            svg.clear();
        },
    
        render : function(){
            var instructs = this.phylogramObject.getDrawInstructs();
            for (var i = 0; i < instructs.length; i++) {
               if(instructs[i].type == 'line'){
                    var line = this.svg.path(["M", instructs[i].x1, instructs[i].y1, "L", instructs[i].x2, instructs[i].y2]).attr(Smits.PhyloCanvas.Render.Style.line);
                } else if(instructs[i].type == 'path'){
                    var path = this.svg.path(instructs[i].path).attr(instructs[i].attr);            
                } else if(instructs[i].type == 'circle'){
                    var path = this.svg.circle(instructs[i].x, instructs[i].y, instructs[i].radius).attr(Smits.PhyloCanvas.Render.Style.text);
                } else {
                    var text = this.svg.text(instructs[i].x, instructs[i].y, instructs[i].text).attr(Smits.PhyloCanvas.Render.Style.text);
                    if(instructs[i].attr){
                        text.attr(instructs[i].attr);
                    }
                    if(instructs[i].rotate){
                        text.rotate(instructs[i].rotate);
                    }
                    
                    var bbox = text.getBBox();
                    var hyp = Math.sqrt( (bbox.height * bbox.height) + (bbox.width * bbox.width) ); // get hypotenuse
                    
                } 
            }
        },
        
        draw : function(instruct){
            var obj, 
                param;

           if(instruct.type == 'line'){
                obj = this.svg.path(["M", instruct.x1, instruct.y1, "L", instruct.x2, instruct.y2]).attr(Smits.PhyloCanvas.Render.Style.line);
            } else if(instruct.type == 'path'){
                obj = this.svg.path(instruct.path).attr(instruct.attr);         
            } else if(instruct.type == 'circle'){
                obj = this.svg.circle(instruct.x, instruct.y, instruct.radius).attr(instruct.attr);
            } else if(instruct.type == 'text'){
                obj = this.svg.text(instruct.x, instruct.y, instruct.text).attr(Smits.PhyloCanvas.Render.Style.text);
                if(instruct.attr){
                    obj.attr(instruct.attr);
                }
                if(instruct.rotate){
                    obj.rotate(instruct.rotate);
                }
                
                var bbox = obj.getBBox();
                param = Math.sqrt( (bbox.height * bbox.height) + (bbox.width * bbox.width) );   // get hypotenuse
            } 

            return [obj, param];
        }   
        
    }
    
};
Smits.PhyloCanvas.Render.Phylogram = function(sSvg, dataObject, refreshFn){

    var svg,
    sParams = Smits.PhyloCanvas.Render.Parameters.Rectangular,  // Easy Reference
    canvasX, canvasY,
    scaleX, scaleY, maxBranch,
    minHeightBetweenLeaves,
    firstBranch = true,
    absoluteY = 0, maxLabelLength = 0,
    outerX, outerY, outerRadius,
    x1, x2, y1, y2, 
    positionX, positionY,
    bufferX, paddingX, paddingY, labelsHold = [],
    
    textPadding = function (y){
        return y + Math.round(y / 4);
    },
    
    rectLinePathArray = function (x1, y1, x2, y2){
        return ["M", x1, y1, "L", x2, y1, "L", x2, y2, "L", x1, y2, "Z"];
    },
    
    recursiveCalculateNodePositions = function (node, positionX){
        if(node.len && firstBranch == false && node.getCountAllChildren() == 0){ 
            absoluteY = Smits.Common.roundFloat(absoluteY + scaleY, 4);
        }
        
        if( node.getCountAllChildren() ){
            var nodeCoords = [], x1,x2,y1,y2;
            if(node.len){ // draw stem
                x1 = positionX;
                x2 = positionX = Smits.Common.roundFloat(positionX + (scaleX * node.len), 4);
                y1 = absoluteY + (node.getMidbranchPosition(firstBranch) * scaleY);
                y2 = y1;
                svg.draw(new Smits.PhyloCanvas.Render.Line(x1, x2, y1, y2));
            }
            
            
            if(node.name){ // draw bootstrap values
                var attr = {};
                attr = Smits.PhyloCanvas.Render.Style.getStyle('bootstrap', 'text');
                if(node.uri) { attr.href = node.uri };
                if(node.description) {attr.title = node.description };
                if(node.level == 0){ 
                    var innerY2 = absoluteY + (node.getMidbranchPosition(firstBranch) * scaleY);
                } else {
                    var innerY2 = y2;
                }
                
                svg.draw(
                    new Smits.PhyloCanvas.Render.Text(
                        (x2 || positionX) + 5, innerY2,
                        node.name,
                        {
                            attr: attr
                        }
                    )
                );          
            }
            
            if( node.getCountAllChildren() ){
                var offspring = node.getRenderChildren();
                for(var i = 0; i < offspring.length; i++){
                    var child = offspring[i];
                    nodeCoords.push(recursiveCalculateNodePositions(child, positionX));
                }
            }
            
            var flatNodeCoords = []; // establish vertical bounds
            for ( var i = 0; i < nodeCoords.length; i++ ){
                if(nodeCoords[i][0]) flatNodeCoords.push(nodeCoords[i][0]);
                if(nodeCoords[i][1]) flatNodeCoords.push(nodeCoords[i][1]);
            }
            var verticalY1 = Math.min.apply(null, flatNodeCoords );
            var verticalY2 = Math.max.apply(null, flatNodeCoords);
            
            // draw vertical
            // hack: little elbows at ends in order to prevent stair-effects at edges
            svg.draw( 
                new Smits.PhyloCanvas.Render.Path( 
                    [
                        "M", positionX + 0.0001, verticalY1,
                        "L", positionX, verticalY1,
                        "L", positionX, verticalY2,
                        "L", positionX + 0.0001, verticalY2
                    ],
                    { attr : Smits.PhyloCanvas.Render.Style.line }
                )               
            );
            
            /* Node Menu */
            var draw = svg.draw(
                new Smits.PhyloCanvas.Render.Circle(positionX, y1, 10,
                    { attr : {fill: 'blue', opacity: 0.001, cursor: 'pointer' } }
                )
            );

            Smits.Common.addRaphEventHandler(
                draw[0], 
                'mouseover', 
                function(o){ 
                    o.raphEl.animate({ opacity: 0.5 }, 200);
                },
                { svg: svg, node: node, x: x2, y: y2, raphEl: draw[0] }
            );          
            Smits.Common.addRaphEventHandler(
                draw[0], 
                'mouseout', 
                function(o){ 
                    o.raphEl.animate({ opacity: 0.001 }, 200);
                },
                { svg: svg, node: node, x: x2, y: y2, raphEl: draw[0] }
            );          
            Smits.Common.addRaphEventHandler(
                draw[0], 
                'click', 
                function(o){ 
                    for( var i = 0; i < o.node.children.length; i++){
                        o.node.children[i].collapsed = true;
                    }
                    o.node.collapsedCap = true;
                    
                    dataObject.getRoot().clearCache();
                    refreshFn();
                },
                { svg: svg, node: node, x: x2, y: y2, raphEl: draw[0] }
            );          

            if (node.color) {
               color = "rgb(" + [ node.color.red[0], node.color.green[0], node.color.blue[0] ].join(",") + ")"
               svg.draw( new Smits.PhyloCanvas.Render.Circle(
                        x1, y1, 15, { attr: { fill: color, opacity: 0.8, cursor: 'hand' } }
                ));
            }
        } else {
            // label
            x1 = positionX;
            x2 = Smits.Common.roundFloat(positionX + (scaleX * node.len), 2);
            y1 = absoluteY;
            y2 = absoluteY;
                
            // preserve for later processing
            node.y = absoluteY;
            labelsHold.push(node);              
                
            svg.draw(new Smits.PhyloCanvas.Render.Line(x1, x2, y1, y2));
            if( sParams.alignRight ){
                svg.draw(
                    new Smits.PhyloCanvas.Render.Path(
                        ["M", x2, y1, "L", sParams.alignPadding + maxBranch, y2],
                        { attr : Smits.PhyloCanvas.Render.Style.connectedDash }
                    )
                );          
            }
            
            if(node.name){
                var attr = {};
                if(node.style){
                    attr = Smits.PhyloCanvas.Render.Style.getStyle(node.style, 'text');
                }
                attr["text-anchor"] = 'start';
                if(node.uri) { attr.href = node.uri };
                if(node.description) {attr.title = node.description };
                
                var draw = svg.draw(
                    new Smits.PhyloCanvas.Render.Text(
                        sParams.alignRight ? maxBranch + sParams.bufferInnerLabels + sParams.alignPadding : x2 + sParams.bufferInnerLabels, y2,
                        node.name,
                        {
                            attr: attr
                        }
                    )
                );              
                maxLabelLength = Math.max(draw[1], maxLabelLength);
                

                // Rollover, Rollout and Click Events
                if(Smits.PhyloCanvas.Render.Parameters.mouseRollOver){
                    Smits.Common.addRaphEventHandler(
                        draw[0], 
                        'mouseover', 
                        Smits.PhyloCanvas.Render.Parameters.mouseRollOver, 
                        { svg: svg, node: node, x: x2, y: y2, textEl: draw[0] }
                    );
                }
                if(Smits.PhyloCanvas.Render.Parameters.mouseRollOut){
                    Smits.Common.addRaphEventHandler(
                        draw[0], 
                        'mouseout', 
                        Smits.PhyloCanvas.Render.Parameters.mouseRollOut, 
                        { svg: svg, node: node, x: x2, y: y2, textEl: draw[0] }
                    );              
                }
                if(Smits.PhyloCanvas.Render.Parameters.onClickAction){
                    Smits.Common.addRaphEventHandler(
                        draw[0], 
                        'click', 
                        Smits.PhyloCanvas.Render.Parameters.onClickAction, 
                        { svg: svg, node: node, x: x2, y: y2, textEl: draw[0] }
                    );              
                }
            }
            
            if( node.collapsedCap ){
                var draw = svg.draw(
                    new Smits.PhyloCanvas.Render.Path(
                        [
                            "M", x2, y1, 
                            "L", x2 + (node.getAllChildrenCombinedMaxLength() * scaleX), y2 + (scaleY / 3),
                            "L", x2 + (node.getImmediateChildrenCombinedMinLength() * scaleX), y2 - (scaleY / 3),
                            x2, y1
                        ],
                        { attr : { stroke: '#000', cursor: 'pointer', fill: '#CCC' } }
                    )
                );          
                Smits.Common.addRaphEventHandler(
                    draw[0], 
                    'click', 
                    function(o){ 
                        for( var i = 0; i < o.node.children.length; i++){
                            o.node.children[i].collapsed = false;
                        }
                        o.node.collapsedCap = false;
                        
                        dataObject.getRoot().clearCache();
                        refreshFn();                
                    },
                    { svg: svg, node: node, x: x2, y: y2, raphEl: draw[0] }
                );  
            }
            
            if(firstBranch){
                firstBranch = false;
            }
        
        }
        
        return [y1, y2];

    },
    
    drawScaleBar = function (){
        y = absoluteY + scaleY;
        x1 = 0;
        x2 = sParams.showScaleBar * scaleX;
        svg.draw(new Smits.PhyloCanvas.Render.Line(x1, x2, y, y));
        svg.draw(new Smits.PhyloCanvas.Render.Text(
            (x1+x2)/2, 
            y-8, 
            sParams.showScaleBar)
        );
    },
    
    renderBinaryChart = function(x, groupName, params){
        var bufferInner = (params && params.bufferInner ? params.bufferInner : 0) | Smits.PhyloCanvas.Render.Parameters.binaryChartBufferInner,
            bufferSiblings = (params && params.bufferSiblings ? params.bufferSiblings * scaleY : 0) | (Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings < 1 ? scaleY * Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings : Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings),     
            thickness = (params && params.thickness ? params.thickness : 0) | Smits.PhyloCanvas.Render.Parameters.binaryChartThickness,
            beginY;
            
        for(var i = 0; i < labelsHold.length; i++){
            var node = labelsHold[i];
            svg.draw(
                new Smits.PhyloCanvas.Render.Path(
                    rectLinePathArray(
                        x + bufferInner,
                        node.y - (scaleY/2) + (bufferSiblings/2),
                        x + bufferInner + thickness, 
                        node.y + (scaleY/2) - (bufferSiblings/2)
                    ),
                    { attr: Smits.PhyloCanvas.Render.Style.getStyle(node.chart[groupName], 'textSecantBg') }
                )
            );          
        }
        return x + bufferInner + thickness;
    },
    
    renderBarChart = function(x, groupName, params){
        var allValues = [], maxValue,
            bufferInner = params && params.bufferInner ? params.bufferInner : 0 | Smits.PhyloCanvas.Render.Parameters.barChartBufferInner,
            height = params && params.height ? params.height : 0 | Smits.PhyloCanvas.Render.Parameters.barChartHeight,
            width = params && params.width ? (params.width < 1 ? scaleY * params.width : params.width ) : 0 | (Smits.PhyloCanvas.Render.Parameters.barChartWidth < 1 ? scaleY * Smits.PhyloCanvas.Render.Parameters.barChartWidth : Smits.PhyloCanvas.Render.Parameters.barChartWidth),
            scaleHeight = 0;
        
        // Need to get max value
        for(var i = 0; i < labelsHold.length; i++){
            allValues.push(labelsHold[i].chart[groupName]);
        }
        maxValue = Math.max.apply(null, allValues);
        scaleHeight = Smits.Common.roundFloat(height / maxValue, 4);
        
        for(var i = 0; i < labelsHold.length; i++){
            var node = labelsHold[i];
            svg.draw(
                    new Smits.PhyloCanvas.Render.Path(
                        rectLinePathArray(
                            x + bufferInner,
                            node.y - (width/2),
                            x + bufferInner + (scaleHeight * node.chart[groupName]), 
                            node.y + (width/2)
                        ),                  
                        { attr: Smits.PhyloCanvas.Render.Style.getStyle(node.chart[groupName], 'barChart') }
                    )
                );                  
        }
        
        return x + bufferInner + height;
    };
    
    
    /* CONSTRUCTOR */
    if(dataObject.getValidate()){   // Validate
        svg.draw(0,0, dataObject.getValidate());
    }

    svg = sSvg;
    var node = dataObject.getRoot();
    var mNewickLen = dataObject.getNewickLen();
    
    canvasX = svg.canvasSize[0];            // Full Canvas Width
    canvasY = svg.canvasSize[1];            // Full Canvas Height
    
    bufferX = sParams.bufferX;
    paddingX = sParams.paddingX;
    paddingY = sParams.paddingY;
    minHeightBetweenLeaves = sParams.minHeightBetweenLeaves;

    absoluteY = paddingY;
    
    scaleX = Math.round((canvasX - bufferX - paddingX*2) / mNewickLen);
    scaleY = Math.round((canvasY - paddingY*2) / (sParams.showScaleBar ? node.getCountAllChildren() : node.getCountAllChildren() - 1 ) );
    if(scaleY < minHeightBetweenLeaves){
        scaleY = minHeightBetweenLeaves;
    }
    maxBranch = Math.round( canvasX - bufferX - paddingX*2 );   
    
    if( (dataObject.binaryCharts && dataObject.binaryCharts.length) || (dataObject.barCharts && dataObject.barCharts.length) ){
        sParams.alignRight = true;
    }
    
    recursiveCalculateNodePositions(node, paddingX);
    
    // Draw Scale Bar
    if(sParams.showScaleBar){
        drawScaleBar();
    }
    
    outerX = maxBranch + maxLabelLength + sParams.bufferInnerLabels;
    // Draw secant highlights
    if(dataObject.binaryCharts && dataObject.binaryCharts.length){
        for( var i = 0; i < dataObject.binaryCharts.length; i++ ){
            outerX = renderBinaryChart(outerX, dataObject.binaryCharts[i].chart, dataObject.binaryCharts[i]);
        }
    }       
    
    // Draw Bar Chart
    if(dataObject.barCharts && dataObject.barCharts.length){
        for( var i = 0; i < dataObject.barCharts.length; i++ ){
            outerRadius = renderBarChart(outerX, dataObject.barCharts[i].chart, dataObject.barCharts[i]);
        }
    }   
    return {    
        getCanvasSize : function(){
            return [canvasX, canvasY];
        },
        getRoot : function(){
            return dataObject.getRoot();
        }

    }
};
/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
var XMLObjectifier = (function() {
    var _clone = function(obj){
        if(!!obj && typeof(obj)==="object"){
            function F(){}
            F.prototype = obj;
            return new F();
        }       
    };
    //Is Numeric check
    var isNumeric = function(s) {
        var testStr = "";
        if(!!s && typeof(s) === "string") { testStr = s; }
        var pattern = /^((-)?([0-9]*)((\.{0,1})([0-9]+))?$)/;
        return pattern.test(testStr);
    };
    var _self = {
    xmlToJSON: function(xdoc) {
        try {
            if(!xdoc){ return null; }
            var tmpObj = {};
                tmpObj.typeOf = "JSXBObject";
            var xroot = (xdoc.nodeType == 9)?xdoc.documentElement:xdoc;
                tmpObj.RootName = xroot.nodeName || "";
            if(xdoc.nodeType == 3 || xdoc.nodeType == 4) {
                return xdoc.nodeValue;
            }
            //Trim function
            function trim(s) {
                return s.replace(/^\s+|\s+$/gm,'');
            }                       
            //Alters attribute and collection names to comply with JS
            function formatName(name) {
                var regEx = /-/g;
                var tName = String(name).replace(regEx,"_");
                return tName;
            }
            //Set Attributes of an object
            function setAttributes(obj, node) {
                if(node.attributes.length > 0) {
                    var a = node.attributes.length-1;
                    var attName;
                    obj._attributes = [];
                    do { //Order is irrelevant (speed-up)
                        attName = String(formatName(node.attributes[a].name));
                        obj._attributes.push(attName);              
                        obj[attName] = trim(node.attributes[a].value);
                    } while(a--);
                }
            }
            
            //Node Prototype
            var _node = (function() {
                    var _self = {
                        activate: function() {
                            var nodes = [];
                            if(!!nodes) {
                                    nodes.getNodesByAttribute = function(attr, obj) {
                                        if(!!nodes && nodes.length > 0) {
                                            var out = [];
                                            var cNode;
                                            var maxLen = nodes.length -1;
                                            try {
                                                do {
                                                    cNode = nodes[maxLen];
                                                    if(cNode[attr] === obj) {
                                                        out.push(cNode);
                                                    }
                                                } while(maxLen--);
                                                out.reverse();
                                                return out;
                                            } catch(e) {return null;}
                                            return null;
                                        }
                                    };
                                    nodes.getNodeByAttribute = function(attr, obj) {
                                        if(!!nodes && nodes.length > 0) {
                                            var cNode;
                                            var maxLen = nodes.length -1;
                                            try {
                                                do {
                                                    cNode = nodes[maxLen];
                                                    if(cNode[attr] === obj) {
                                                        return cNode;
                                                    }
                                                } while(maxLen--);
                                            } catch(e) {return null;}
                                            return null;
                                        }
                                    };
                                    nodes.getNodesByValue = function(obj) {
                                        if(!!nodes && nodes.length > 0) {
                                            var out = [];
                                            var cNode;
                                            var maxLen = nodes.length -1;
                                            try {
                                                do {
                                                    cNode = nodes[maxLen];
                                                    if(!!cNode.Text && cNode.Text === obj) {
                                                        out.push(cNode);
                                                    }
                                                } while(maxLen--);
                                                return out;
                                            } catch(e) {return null;}
                                            return null;
                                        }
                                    };
                                    nodes.contains = function(attr, obj) {
                                        if(!!nodes && nodes.length > 0) {
                                            var maxLen = nodes.length -1;
                                            try {
                                                do {
                                                    if(nodes[maxLen][attr] === obj) {
                                                        return true;
                                                    }
                                                } while(maxLen--);
                                            } catch(e) {return false;}
                                            return false;
                                        }
                                    };
                                    nodes.indexOf = function(attr, obj) {
                                        var pos = -1;
                                        if(!!nodes && nodes.length > 0) {
                                            var maxLen = nodes.length -1;
                                            try {
                                                do {
                                                    if(nodes[maxLen][attr] === obj) {
                                                        pos = maxLen;
                                                    }
                                                } while(maxLen--);
                                            } catch(e) {return -1;}
                                            return pos;
                                        }
                                    };
                                    nodes.SortByAttribute = function(col, dir) {
                                        if(!!nodes && nodes.length > 0) {               
                                            function getValue(pair, idx) {
                                                var out = pair[idx];
                                                out = (bam.validation.isNumeric(out))?parseFloat(out):out;
                                                return out;
                                            }
                                            function sortFn(a, b) {
                                                var tA, tB;
                                                tA = getValue(a, col);
                                                tB = getValue(b, col);
                                                var res = (tA<tB)?-1:(tB<tA)?1:0;
                                                if(!!dir) {
                                                    res = (dir.toUpperCase() === "DESC")?(0 - res):res;
                                                }
                                                return res;
                                            }
                                            nodes.sort(sortFn);
                                        }
                                    };
                                    nodes.SortByValue = function(dir) {
                                        if(!!nodes && nodes.length > 0) {
                                            function getValue(pair) {
                                                var out = pair.Text;
                                                out = (bam.validation.isNumeric(out))?parseFloat(out):out;
                                                return out;
                                            }
                                            function sortFn(a, b) {
                                                var tA, tB;
                                                tA = getValue(a);
                                                tB = getValue(b);
                                                var res = (tA<tB)?-1:(tB<tA)?1:0;
                                                if(!!dir) {
                                                    res = (dir.toUpperCase() === "DESC")?(0 - res):res;
                                                }
                                                return res;
                                            }
                                            nodes.sort(sortFn);
                                        }
                                    };
                                    nodes.SortByNode = function(node, dir) {
                                        if(!!nodes && nodes.length > 0) {
                                            function getValue(pair, node) {
                                                var out = pair[node][0].Text;
                                                out = (bam.validation.isNumeric(out))?parseFloat(out):out;
                                                return out;
                                            }
                                            function sortFn(a, b) {                                     
                                                var tA, tB;
                                                tA = getValue(a, node);
                                                tB = getValue(b, node);
                                                var res = (tA<tB)?-1:(tB<tA)?1:0;
                                                if(!!dir) {
                                                    res = (dir.toUpperCase() === "DESC")?(0 - res):res;
                                                }
                                                return res;
                                            }
                                            nodes.sort(sortFn);
                                        }
                                  };
                            }
                            return nodes;
                        }
                    };
                    return _self;
            })();
            //Makes a new node of type _node;
            var makeNode = function() {
                var _fn = _clone(_node);                    
                return _fn.activate();
            }
            //Set collections
            function setHelpers(grpObj) {
                //Selects a node withing array where attribute = value
                grpObj.getNodeByAttribute = function(attr, obj) {
                    if(this.length > 0) {
                        var cNode;
                        var maxLen = this.length -1;
                        try {
                            do {
                                cNode = this[maxLen];
                                if(cNode[attr] == obj) {
                                    return cNode;
                                }
                            } while(maxLen--);
                        } catch(e) {return false;}
                        return false;
                    }
                };
                
                grpObj.contains = function(attr, obj) {
                    if(this.length > 0) {
                        var maxLen = this.length -1;
                        try {
                            do {
                                if(this[maxLen][attr] == obj) {
                                    return true;
                                }
                            } while(maxLen--);
                        } catch(e) {return false;}
                        return false;
                    }
                };
                
                grpObj.indexOf = function(attr, obj) {
                    var pos = -1;
                    if(this.length > 0) {
                        var maxLen = this.length -1;
                        try {
                            do {
                                if(this[maxLen][attr] == obj) {
                                    pos = maxLen;
                                }
                            } while(maxLen--);
                        } catch(e) {return -1;}
                        return pos;
                    }
                };
                
                grpObj.SortByAttribute = function(col, dir) {
                    if(this.length) {               
                        function getValue(pair, idx) {
                            var out = pair[idx];
                            out = (isNumeric(out))?parseFloat(out):out;
                            return out;
                        }
                        function sortFn(a, b) {
                            var res = 0;
                            var tA, tB;                     
                            tA = getValue(a, col);
                            tB = getValue(b, col);
                            if(tA < tB) { res = -1; } else if(tB < tA) { res = 1; }
                            if(dir) {
                                res = (dir.toUpperCase() == "DESC")?(0 - res):res;
                            }
                            return res;
                        }
                        this.sort(sortFn);
                    }
                };
                
                grpObj.SortByValue = function(dir) {
                    if(this.length) {
                        function getValue(pair) {
                            var out = pair.Text;
                            out = (isNumeric(out))?parseFloat(out):out;
                            return out;
                        }
                        function sortFn(a, b) {
                            var res = 0;
                            var tA, tB;
                            tA = getValue(a);
                            tB = getValue(b);
                            if(tA < tB) { res = -1; } else if(tB < tA) { res = 1; }
                            if(dir) {
                                res = (dir.toUpperCase() == "DESC")?(0 - res):res;
                            }
                            return res;
                        }
                        this.sort(sortFn);
                    }
                };
                
                grpObj.SortByNode = function(node, dir) {
                    if(this.length) {
                        function getValue(pair, node) {
                            var out = pair[node][0].Text;
                            out = (isNumeric(out))?parseFloat(out):out;
                            return out;
                        }
                        function sortFn(a, b) {
                            var res = 0;
                            var tA, tB;
                            tA = getValue(a, node);
                            tB = getValue(b, node);
                            if(tA < tB) { res = -1; } else if(tB < tA) { res = 1; }
                            if(dir) {
                                res = (dir.toUpperCase() == "DESC")?(0 - res):res;
                            }
                            return res;
                        }
                        this.sort(sortFn);
                    }
                };
            }
            //Recursive JSON Assembler
            //Set Object Nodes
            function setObjects(obj, node) {
                var elemName;   //Element name
                var cnode;  //Current Node
                var tObj;   //New subnode
                var cName = "";
                if(!node) { return null; }              
                //Set node attributes if any
                if(node.attributes.length > 0){setAttributes(obj, node);}               
                obj.Text = "";
                if(node.hasChildNodes()) {
                    var nodeCount = node.childNodes.length - 1; 
                    var n = 0;
                    do { //Order is irrelevant (speed-up)
                        cnode = node.childNodes[n];
                        switch(cnode.nodeType) {
                            case 1: //Node
                            //Process child nodes
                            obj._children = [];
                            //SOAP XML FIX to remove namespaces (i.e. soapenv:)
                            elemName = (cnode.localName)?cnode.localName:cnode.baseName;
                            elemName = formatName(elemName);
                            if(cName != elemName) { obj._children.push(elemName); }
                                //Create sub elemns array
                                if(!obj[elemName]) {
                                    obj[elemName] = []; //Create Collection
                                }
                                tObj = {};
                                obj[elemName].push(tObj);
                                if(cnode.attributes.length > 0) {
                                    setAttributes(tObj, cnode);
                                }
                                //Set Helper functions (contains, indexOf, sort, etc);
                                if(!obj[elemName].contains) {
                                    setHelpers(obj[elemName]);
                                }   
                            cName = elemName;
                            if(cnode.hasChildNodes()) {
                                setObjects(tObj, cnode); //Recursive Call
                            }
                            break;
                            case 3: //Text Value
                            obj.Text += trim(cnode.nodeValue);
                            break;
                            case 4: //CDATA
                            obj.Text += (cnode.text)?trim(cnode.text):trim(cnode.nodeValue);
                            break;
                        }
                    } while(n++ < nodeCount);
                }
            }           
            //RUN
            setObjects(tmpObj, xroot);
            //Clean-up memmory
            xdoc = null;
            xroot = null;
            return tmpObj;  
        } catch(e) {
                return null;    
        }   
    },

    //Converts Text to XML DOM
    textToXML: function(strXML) {
        var xmlDoc = null;
        try {
            xmlDoc = (document.all)?new ActiveXObject("Microsoft.XMLDOM"):new DOMParser();
            xmlDoc.async = false;
        } catch(e) {throw new Error("XML Parser could not be instantiated");}
        var out;
        try {
            if(document.all) {
                out = (xmlDoc.loadXML(strXML))?xmlDoc:false;
            } else {        
                out = xmlDoc.parseFromString(strXML, "text/xml");
            }
        } catch(e) { throw new Error("Error parsing XML string"); }
        return out;
    }
    };
    return _self;
})();

