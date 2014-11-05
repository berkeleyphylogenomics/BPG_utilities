/* Contains functions dealing with the arborjs setup of the graph in phylo4j explorer. */
// for the raphael paper
var graphPaper;

// For the arbor js rendering functions
var explorerRenderer = function(canvas) {
    // canvas padding
    var padding = 0;
    // get height, width
    var h = $j('#'+canvas).height(), w = $j('#'+canvas).width();
    // make a new raphael canvas
    graphPaper = Raphael(canvas, w - padding, h - padding);
    // arbor js particle system
    var particleSystem;
    // radius of the node
    var nodeRadius = 10;
    // number of pixels between bezier control points
    var bezierGrade = 10;
    // height of the graph
    var graphHeight = 450;

    var that = {
        init: function(system) {
            particleSystem = system;
            particleSystem.screenSize(w,h);
            particleSystem.screenPadding(padding);
            // bind the resize handlers
            //$j(' .explorer-display ').resize(that.resize);
            $j(window).resize(that.resize);
            that.resize();
            that.initMouseHandling();
        },
        redraw: function() {
            // clear the paper
            graphPaper.clear();
            // does a redraw of the svg paper
            var edgesDone = [];
            // this function gives an iterator for each node
            particleSystem.eachNode(function(node, pt) {
                var thisNode = graphPaper.circle(pt.x,pt.y,nodeRadius).attr(node.data.attr);
            });
            // iterates over every edge
            particleSystem.eachEdge(function(edge,pt1,pt2) {
                // have we already done this edge?
                if($j.inArray(edge, edgesDone) == -1) {
                    var theseEdges = particleSystem.getEdges(edge.source,edge.target).concat(particleSystem.getEdges(edge.target,edge.source));
                    var thisBezierHeightArray = returnBezierHeightArray(theseEdges.length);
                    var conSource = particleSystem.toScreen(edge.target.p), conTarget = particleSystem.toScreen(edge.source.p);
                    var edgePoints = circleLineIntersectionPoints(conSource, conTarget, nodeRadius);
                    $j.each(theseEdges, function(index,edge) {
                        var bezierPoint = bezierControlPoint(thisBezierHeightArray[index]*bezierGrade,
                                                conSource, conTarget);
                        // change the path
                        var thisEdge = graphPaper.path('M'+edgePoints.p1.x.toString()+','+
                                                edgePoints.p1.y.toString()+' Q'+
                                                bezierPoint.x.toString()+','+
                                                bezierPoint.y.toString()+' '+
                                                edgePoints.p2.x.toString()+','+
                                                edgePoints.p2.y.toString()).attr(edge.data.attr); 
                        edgesDone.push(edge);
                    });
                }
            }); 
        },
        resize: function() {
            w = $j('#'+canvas).width();
            //h = $j('#'+canvas).height();
            h = graphHeight;
            graphPaper.setSize(w,h);
            particleSystem.screenSize(w,h);
            that.redraw();
        },
        initMouseHandling: function() {
            
        }
    };
    return that;
}
var particleSystemProperties = {
    // properties for the forces in the force-directed graph
    repulsion: 1,
    stiffness: 2,
    friction: 1,
    gravity: true,
    fps:55,
    dt:0.02,
    precision:0.6
};
// 
var pSys;

function returnBezierHeightArray(n) {
    // takes an integer n (number of edges between 2 nodes), and returns
    // an array of length n for the bezier control point heights
    var ret = [];
    var range = Math.floor(n / 2);
    for(var i=-1*range;i<=range;i++) {
        if(!((i == 0)&&(n % 2 == 0))) {ret.push(i);}
    }
    return ret;
}

function bezierControlPoint(h,sp,tp) {
    // given a height above the line, return the coordinates of a point perpendicular
    // to the line between the nodes
    if(h == 0) {return arbor.Point(Math.round((sp.x+tp.x)/2), Math.round((sp.y+tp.y)/2));}
    // address infinite slope condition
    if(sp.y == tp.y) {return arbor.Point(Math.round((sp.x+tp.x)/2), sp.y+h);}
    var thisH = Math.abs(h);
    var edgeSlope = (-1.0 * (tp.x - sp.x) / (tp.y - sp.y));
    var xm = Math.round((sp.x+tp.x)/2), ym = Math.round((sp.y+tp.y)/2);
    var xc;
    if(h < 0) {xc = 1.0*xm - Math.sqrt(1.0*Math.pow(thisH,2)/(1 + Math.pow(edgeSlope,2)));}
    else {xc = 1.0*xm + Math.sqrt(1.0*Math.pow(thisH,2)/(1 + Math.pow(edgeSlope,2)));}
    var yc = edgeSlope * (xc - xm) + ym;
    return arbor.Point(Math.round(xc),Math.round(yc));
}

function circleLineIntersectionPoints(pt1,pt2,r) {
    // returns the two points on the node edge where the edge should actually begin
    // given the 2 center points.
    var edgeSlope = (1.0 * (pt2.y - pt1.y)) / (pt2.x -pt1.x);
    var magicFactorX = Math.sqrt(Math.pow(r, 2) / (1 + Math.pow(edgeSlope, 2)));
    var magicFactorY = edgeSlope * magicFactorX;
    var pt1x_a = Math.round(magicFactorX + pt1.x), pt1x_b = Math.round(pt1.x - magicFactorX);
    var pt2x_a = Math.round(magicFactorX + pt2.x), pt2x_b = Math.round(pt2.x - magicFactorX);
    var pt1y_a = Math.round(magicFactorY + pt1.y), pt1y_b = Math.round(pt1.y - magicFactorY);
    var pt2y_a = Math.round(magicFactorY + pt2.y), pt2y_b = Math.round(pt2.y - magicFactorY);
    var chosenX1, chosenX2, chosenY1, chosenY2;
    if(pt1.x<=pt2.x) {
        chosenX1 = pt1x_a;
        chosenX2 = pt2x_b;
    }
    else {
        chosenX1 = pt1x_b;
        chosenX2 = pt2x_a;
    }
    if(pt1.y<=pt2.y) {
        chosenY1 = pt1y_a;
        chosenY2 = pt2y_b;
    }
    else {
        chosenY1 = pt1y_b;
        chosenY2 = pt2y_a;
    }
    return {p1:arbor.Point(chosenX1,chosenY1),p2:arbor.Point(chosenX2,chosenY2)}
}

function updateExplorerGraphDisplay() {
    // this function will update the display of the graph according to the new options
    // maybe turn on loading indicators, disable buttons?
    //pSys.stop();

    // first get rid of the previous system
    pSys.eachEdge(function(edge, pt1, pt2) {pSys.pruneEdge(edge);});
    pSys.eachNode(function(node,pt) {pSys.pruneNode(node);});

    graphPaper.remove();
    //graphPaper = ''
    // pSys = '';
    //$j(' .explorer-graph-display canvas ').remove();

    // Get the graph topology data from the server
    $j.ajax('/api/phylo4j/explorer/' + nodeID.toString() + '/graph_data/', {
        data: {
            'depth': $j(' #graph-depth-control ').val(),
            'relationship_types': getRelationshipTypes()
        },
        success: function(data) {
            thisGraphData = data;
            if ('error' in data) {$j(' .explorer-graph-display ').html(data['error']);}
            else {
                // add the canvas back
                //$j(' .explorer-graph-display ').html("<canvas></canvas>");
                // load the graph system
                pSys = arbor.ParticleSystem(particleSystemProperties);
                //pSys.parameters({gravity:true});
                pSys.renderer = explorerRenderer('explorer-graph-display');
                pSys.graft(data);
                pSys.start();
            }
        }
    });
}

