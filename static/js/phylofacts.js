// Eventually, we should have a single global variable (say, BPG) and 
// all functions should be attributes of that global variable.

// This function makes several input boxes on index.html visible.
var show_input_boxes = function(){
    jQuery("#controlPanel ul li a").removeClass("light");
    jQuery("#" + this.id).addClass('light');
    jQuery("#menu-options").hide();
    jQuery(".panel").hide();
    jQuery("#" + this.id + "-panel").show();
    jQuery(".input-box").focus();
    return false;
}
  
var EXAMPLE_VALUES = {"phylofacts-family-search" : "bpg0240116",
    "sequence-search" : ">sp|O14727|APAF_HUMAN Apoptotic protease-activating factor 1 OS=Homo sapiens GN=APAF1 PE=1 SV=2\nMDAKARNCLLQHREALEKDIKTSYIMDHMISDGFLTISEEEKVRNEPTQQQRAAMLIKMILKKDNDSYVSFYNALLHEGYKDLAALLHDGIPVVSSSSGKDSVSGITSYVRTVLCEGGVPQRPVVFVTRKKLVNAIQQKLSKLKGEPGWVTIHGMAGCGKSVLAAEAVRDHSLLEGCFPGGVHWVSVGKQDKSGLLMKLQNLCTRLDQDESFSQRLPLNIEEAKDRLRILMLRKHPRSLLILDDVWDSWVLKAFDSQCQILLTTRDKSVTDSVMGPKYVVPVESSLGKEKGLEILSLFVNMKKADLPEQAHSIIKECKGSPLVVSLIGALLRDFPNRWEYYLKQLQNKQFKRIRKSSSYDYEALDEAMSISVEMLREDIKDYYTDLSILQKDVKVPTKVLCILWDMETEEVEDILQEFVNKSLLFCDRNGKSFRYYLHDLQVDFLTEKNCSQLQDLHKKIITQFQRYHQPHTLSPDQEDCMYWYNFLAYHMASAKMHKELCALMFSLDWIKAKTELVGPAHLIHEFVEYRHILDEKDCAVSENFQEFLSLNGHLLGRQPFPNIVQLGLCEPETSEVYQQAKLQAKQEVDNGMLYLEWINKKNITNLSRLVVRPHTDAVYHACFSEDGQRIASCGADKTLQVFKAETGEKLLEIKAHEDEVLCCAFSTDDRFIATCSVDKKVKIWNSMTGELVHTYDEHSEQVNCCHFTNSSHHLLLATGSSDCFLKLWDLNQKECRNTMFGHTNSVNHCRFSPDDKLLASCSADGTLKLWDATSANERKSINVKQFFLNLEDPQEDMEVIVKCCSWSADGARIMVAAKNKIFLFDIHTSGLLGEIHTGHHSTIQYCDFSPQNHLAVVALSQYCVELWNTDSRSKVADCRGHLSWVHGVMFSPDGSSFLTSSDDQTIRLWETKKVCKNSAVMLKQEVDVVFQENEVMVLAVDHIRRLQLINGRTGQIDYLTEAQVSCCCLSPHLQYIAFGDENGAIEILELVNNRIFQSRFQHKKTVWHIQFTADEKTLISSSDDAEIQVWNWQLDKCIFLRGHQETVKDFRLLKNSRLLSWSFDGTVKVWNIITGNKEKDFVCHQGTVLSCDISHDATKFSSTSADKTAKIWSFDLLLPLHELRGHNGCVRCSAFSVDSTLLATGDDNGEIRIWNVSNGELLHLCAPLSEEGAATHGGWVTDLCFSPDGKMLISAGGYIKWWNVVTGESSQTFYTNGTNLKKIHVSPDFKTYVTVDNLGILYILQTLE",
   "orthology-search" : "AAT_ECOLI",
   "accession-search" : "APAF_HUMAN",
   "phylofacts-pfam-search" : "PF00001",
   "phylofacts-biocyc-search" : "3.1.1.1"}

// This function is used to display example queries in input boxes based on the values in EXAMPLE_VALUES. 
var show_example_value = function(obj_id){
    jQuery("#" + obj_id + "-example").click(function(){
        jQuery("#" + obj_id + "-field").val(EXAMPLE_VALUES[obj_id]);
        return false;
    });
}

var validate = function(form_id){
    jQuery("#" + form_id + "-form").submit(function(){
        if (jQuery("#" + form_id + "-form .input-box").val() == ""){
            jQuery("#" + form_id + "-errors").text("This field cannot be blank.").show();
            jQuery("#" + form_id + "-field").focus();
            return false;
        }
    });
};

// -------------------- end of function definitions -----------------------

jQuery(document).ready(function(){
    validate("phylofacts-family-search");
    validate("sequence-search");
    validate("orthology-search");
    validate("accession-search");
    validate("phylofacts-pfam-search");
    validate("phylofacts-biocyc-search");
});

jQuery(document).ready(function(){
    jQuery(".close-panel").click(function(){
        jQuery(".panel").hide();
        jQuery("#menu-options").show();
        return false;
    });
});

// Input boxes on the home page rightPanel are hidden by default. This code attached
// event-handlers to make them visible.
jQuery(document).ready(function(){
    jQuery("#sequence-search").click(show_input_boxes);
    jQuery("#phylofacts-pfam-search").click(show_input_boxes);
    jQuery("#phylofacts-biocyc-search").click(show_input_boxes);
    jQuery("#accession-search").click(show_input_boxes);
    jQuery("#phylofacts-family-search").click(show_input_boxes);
    jQuery("#orthology-search").click(show_input_boxes);
    jQuery("#download").click(show_input_boxes);
});

// This code attaches a click event-handler to each 'Example' button that is found on any of 
// the forms that are dynamically displayed on the home page.
jQuery(document).ready(function(){
    for (var i in EXAMPLE_VALUES){
        show_example_value(i);
    }
});

// This code makes it possible to display the full text of the intro on the index page.
// It also permits the user to revert to a partial display of the intro.
jQuery(document).ready(function(){
    jQuery("#show-full-description-link").click(function(event){
        event.preventDefault();
        jQuery("#hidden-description").show();
        jQuery(this).hide();
    });

    jQuery("#hide-full-description-link").click(function(event){
        event.preventDefault();
        jQuery("#hidden-description").hide();
        jQuery("#show-full-description-link").show();
    });
});


// TODO: refactor this and make it DRY
// This function is useful when there are a lot of GO annotations to 
// display. With the help of 'more' and 'less' anchors, it 
// helps expand or contract the list of annotations displayed.
jQuery(document).ready(function(){
  jQuery(".show-more-go-data").click(function(){
    var process = this.id.split("-")[1]

    jQuery("#more-" + process).hide();
    jQuery(".hidden-go-data-" + process).show();
    jQuery("#less-" + process).show();
    return false;
  });
  
  jQuery(".show-less-go-data").click(function(){
    var process = this.id.split("-")[1]
    jQuery("#less-" + process).hide();
    jQuery(".hidden-go-data-" + process).hide();
    jQuery("#more-" + process).show();
    return false;
  });
});


// This function shows the form on the family page to do another
// search.
jQuery(document).ready(function(){
  jQuery("#another-search").click(function(){
    jQuery(this).hide();
    jQuery("form").show();
    return false;
  });

  jQuery("#hide-form").click(function(){
    jQuery("form").hide();
    jQuery("#another-search").show();
    return false;
  });
});


// This function creates the popup to show the Jalview applet
jQuery(document).ready(function(){
  jQuery("#jalview-applet-for-treenode").click(function() {
    var popup = window.open('', '', "height=200, width=350");  
    if (window.focus) {popup.focus();}
    var treenode = jQuery("#jalview-applet-for-treenode").attr('treenode');
    var tmp = popup.document;
    tmp.write("<html><head>");
    tmp.write("<title>Jalview</title>");
    tmp.write("</head><body>");
    tmp.write("<div id='alignmentinfo'>");
    tmp.write("<applet code='jalview.bin.JalviewLite' width='320' height='35' archive='/static/java/jalviewApplet.jar'>");
    tmp.write("<param name='file' value='/phylofacts/tree_node_view/" + treenode + "/alignment/' />");
    tmp.write("<param name='label' value='Jalview Alignment Viewer' />");
    tmp.write("<param name='userDefinedColourr' value='C=yellow; R,K,H=FF5555; D,E=5555FF' />");
    tmp.write("<param name='showFullId' value='false' />");
    tmp.write("<param name='showTreeDistances' value='false' />");
    tmp.write("<param name='showTreeBootstraps' value='false' />");
    tmp.write("<param name='RGB' value='DDDDDD' />");
    tmp.write("<param name='linkLabel_1' value='SRS' />");
    tmp.write("<param name='linkUrl_1' value='http://srs.ebi.ac.uk/srs7bin/cgi-bin/wgetz?-e+[uniprot-all:jQuerySEQUENCE_IDjQuery]+-vn+2' />");
    tmp.write("<param name='linkLabel_2' value='Uniprot' />");
    tmp.write("<param name='linkUrl_2' value='http://us.expasy.org/cgi-bin/niceprot.pl?jQuerySEQUENCE_IDjQuery' />");
    tmp.write("<param name='APPLICATION_URL' value='http://www.jalview.org/services/launchApp' />");
    tmp.write("The <a href='http://www.jalview.org/' target='_blank'>Jalview</a> applet failed to load; please <a href='http://java.com/en/download/'>download Java</a>.");
    tmp.write("</applet>");
    tmp.write("<p>Powered by Jalview <a href='http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp033' target='_blank'>doi: 10.1093/bioinformatics/btp033</a></p>");
    tmp.write("</div>");
    tmp.write("</body></html>");
    tmp.close();
    return false;
  });
  
  jQuery("#show-lineage").click(function(event){
    event.preventDefault();
    jQuery("#lineage").show();
    jQuery(this).hide();
    jQuery("#hide-lineage").show();
  });

  jQuery("#hide-lineage").click(function(event){
    event.preventDefault();
    jQuery("#lineage").hide();
    jQuery(this).hide();
    jQuery("#show-lineage").show();
  });
});

// You know what, I like jquery.  I was going to write this all in YUI kind of style, but sans YUI,
// but we use jquery, not yui, and jquery has its own shiny widget creation jazz.
//
// Alignment Map
//
// Looks like BLAST's alignment thing
(function ($) {
    $.widget("ui.alignmentmap", {
        // These options will be used as defaults
        options: {
            masterLength: 1,
            // Length of the master sequence
            matches: [
            // Additional fields are allowed, but this is the minimum
            // { alignment_match_from: 10, alignment_match_to: 30, evalue: 1e-10 }
            ],
            rampUpperColor: [0, 104, 204],
            rampLowerColor: [219, 0, 153],
            ramp: function r(evalue) { // color ramp function
                // Everything above this is blue
                var upper_bound = 1e-4;
                // Everything below this is red
                var lower_bound = 1e-30;

                var upper_color = this.rampUpperColor;
                var lower_color = this.rampLowerColor;

                if (evalue >= upper_bound) {
                    return upper_color;
                }
                if (evalue <= lower_bound) {
                    return lower_color;
                }

                return _.map(_.zip(upper_color, lower_color), function (pair) {
                    var dv = pair[0] - pair[1];
                    var v = (Math.log(evalue) - Math.log(lower_bound)) / (Math.log(upper_bound) - Math.log(lower_bound));
                    return pair[1] + v * dv;
                });
            },
            rowHeight: 5,
            rowSpacing: 3,
            matchMouseclick: function (match) {
                return function() {}
            },
            matchMouseover: function (match) {
                return function () {}
            },
            matchMouseup: function (match) {
                return function () {}
            },
            matchMouseout: function (match) {
                return function () {}
            },
            matchCallback: function (match, node) {} // called after the match is rendered
        },

        // Set up the widget
        _create: function () {
            // Should this not use [0]?
            this.paper = Raphael(this.element[0], $(this.element[0]).width(), $(this.element[0]).height());
            this.refresh();
        },

        refresh: function () {
            // If we aren't visible, sod off
            if ($(this.element[0]).is(":hidden")) {
                return
            }
            this.paper.clear();
            var width = $(this.element[0]).width() - 10;
            var height = this.options.matches.length * (this.options.rowHeight + this.options.rowSpacing) + // matches
            this.options.rowHeight * 4; //header and bottom area
            var scale = width / this.options.masterLength;

            this.paper.setSize(width + 10, height);

            // Master Sequence
            this.paper.rect(0, 0, width, this.options.rowHeight * 2);
            this.paper.text(10, this.options.rowHeight, "1");
            this.paper.text(width - 20, this.options.rowHeight, this.options.masterLength);


            _.each(this.options.matches, function (match, index) {
                var match_rect = this.paper.rect(
                match.alignment_match_from * scale, (index + 2) * (this.options.rowHeight + this.options.rowSpacing), (match.alignment_match_to - match.alignment_match_from) * scale, this.options.rowHeight).attr({
                    fill: "rgb(" + this.options.ramp(match.evalue).join(',') + ")",
                    'stroke': 'none', 'cursor':'pointer'
                });
                match_rect.mouseover(this.options.matchMouseover(match));
                match_rect.mouseup(this.options.matchMouseup(match));
                match_rect.mouseout(this.options.matchMouseout(match));
                match_rect.click(this.options.matchMouseclick(match));
                this.options.matchCallback(match, match_rect);
            }, this);
        }
    });
}(jQuery));

//
// PhyloTree
// 
// It's a tree viewer for phyloXML
(function ($) {
    // Maybe this should be monkey patched, maybe not
    function add_tooltip_to_node(r, content) {
        $(r.node).qtip({
            content: content,
            position: {
                my: 'bottom right',
                at: 'top center',
                viewport: $(window),
            },
            hide: {
                fixed: true
            },
            style: {
                classes: 'ui-tooltip-shadow',
                widget: true,
                width: '400px'
            }
        });

    };


    $.widget("ui.phylotree", {
        // These options will be used as defaults
        options: {
            data: {},
            rowHeight: 25,
            showScale: true, // If we want to show the tree length ruler thing
            edgeStyle: {
                stroke: "rgb(90,0,200)"
            },
            nodeStyle: {
                fill: "rgba(0, 104, 204,0.7)"
            },
            // The line from the root to highlighted nodes
            highlightStrokeStyle: {
                "stroke": "rgba(200,0,10,0.5)", 
                "stroke-width": 5
            },
            highlightNodeCallback: function() {},
            highlightPathCallback: function() {},
            pathCallback: function() {},
            // Area to the left, mainly used to take into consideration the line width
            paddingLeft: 2,
	        // Number of characters to limit the description field in the leaves to
            descriptionSubstringLimit: 25,
            leafLabelAttributes: ["name" ],
            labelMakers: {
                // TODO abstract and clean this up.
                description: function (node, x, y) {
                    var sequence_text = '';
                    if (node.getParameters().sequences && node.getParameters().sequences[0] && node.getParameters().sequences[0].name && node.getParameters().sequences[0].name[0]) {
                        sequence_text = node.getParameters().sequences[0].name[0].Text || '';
                    }
                    //// The condition below is for summary PHOGs or other Phylofacts description of summary leaves.
                    else {
                        if (node.getParameters().properties) {
                            var descriptions = _.filter(node.getParameters().properties, function (p) {
                                return p.ref === "PhyloFacts::Description";
                            });
                            if (descriptions.length) {
                                // Only does the first one.  this is stupid
                                sequence_text = descriptions[0].Text;
                            }
                        }
                    }
                    if (sequence_text.length) {
                        var reduced_sequence_text = sequence_text.substring(0, sequence_text.lastIndexOf(' ', this.options.descriptionSubstringLimit)) + ' ...';
                        var r = this.paper.text(x + 2, y, reduced_sequence_text).attr({
                            'text-anchor': 'start'
                        });
                        add_tooltip_to_node(r, sequence_text);
                        return r.getBBox().x2 + 1;
                    } else {
                        return x;
                    }
                },
                ec: function (node, x, y) {
                    var ec = _.filter(node.getParameters().sequences[0].annotation, function (a) {
                        return a.ref.split(":")[0] === "EC";
                    });
                    return _.reduce(ec, function (m, a) {
                        return this.paper.text(m + 2, y, a.ref).attr({
                            'text-anchor': 'start',
                            'fill': '#06C',
                            'href' : "http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=" + a.ref.split(":")[1]
                        }).getBBox().x2
                    }, x, this);
                },
                go_exp: function (node,x, y) {
                    var go = _.filter(node.getParameters().sequences[0].annotation, function (a) {
                        return a.ref.split(":")[0] === "GO";
                    });
                    var go_exp = _.filter(go, function (a) {
                        return _.include(["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"], a.evidence);
                    });
                    if (go_exp.length) {
                        var r = this.paper.image('/static/img/icons/icon_flask_12.png', x + 2, y - 4, 8, 8).attr({
                            'text-anchor': 'start'
                        });
                        var content = _.map(go_exp, function(a) { return  a.desc[0].Text + ' (' + a.evidence + ')'; }).join('<br>');
                        add_tooltip_to_node(r, content);
                        return r.getBBox().x2;
                    } else {
                        return x;
                    }
                },
                go: function (node,x, y) {
                    var go = _.filter(node.getParameters().sequences[0].annotation, function (a) {
                        return a.ref.split(":")[0] === "GO";
                    });
                    if (go.length) {
                        var r = this.paper.image('/static/img/icons/icon_flask_12.png', x + 2, y - 4, 8, 8).attr({
                            'text-anchor': 'start'
                        });
                        var content = _.map(go, function(a) { return  a.desc[0].Text + '(' + a.evidence + ')'; }).join('<br>');
                        add_tooltip_to_node(r, content);
                        return r.getBBox().x2;
                    } else {
                        return x;
                    }
                },
                sfld: function (node,x, y) {
                    var sfld = _.filter(node.getParameters().sequences[0].annotation, function (a) {
                        return a.ref.split(":")[0] === "SFLD";
                    });
                    if (sfld.length) {
                        var r = this.paper.image('/static/img/icons/sfld.png', x + 2, y - 6, 24, 10).attr({
                            'text-anchor': 'start'
                        });
                        console.log(sfld);
                        var content = _.map(sfld, function(a) { return  a.Text; }).join('<br>');
                        add_tooltip_to_node(r, content);
                        return r.getBBox().x2;
                    } else {
                        return x;
                    }
                },
                name: function (node, x, y) {
                    return this.paper.text(x, y, node.getParameters().name).attr({
                         'text-anchor': 'start',
                         'fill': '#06C',
                         'href' : (node.getParameters().name.match(/^PHOG.*/) ? "/phog/" : "/phylofacts/sequence/UniProt/" ) + node.getParameters().name
                    }).getBBox().x2;
                },
                taxonomy: function (node, x, y) {
                    if (node.getParameters().taxonomies && node.getParameters().taxonomies[0] && node.getParameters().taxonomies[0].scientific_name && node.getParameters().taxonomies[0].scientific_name[0]) {
                        // TODO check to make sure ID is there
                        var r = this.paper.text(x + 2, y, node.getParameters().taxonomies[0].scientific_name[0].Text).attr({
                            'text-anchor': 'start',
                            'font-style': 'italic',
                            'fill': '#06C',
                            'href': 'http://www.uniprot.org/taxonomy/' + node.getParameters().taxonomies[0].id[0].Text
                        })
                        var content = node.getParameters().taxonomies[0].scientific_name[0].Text + "<br>";
                        if ( node.getParameters().taxonomies[0].common_name[0].Text ) {
                            content += node.getParameters().taxonomies[0].common_name[0].Text + "<br>";
                        }
                        content += "<i>Click for UniProt Taxonomy Page</i>";
                        add_tooltip_to_node(r, content);
                        return r.getBBox().x2;
                    }
                    return x;
                }
            },
            renderLabel: function (node) {
                var x = this._x(node);
                var y = this._y(node);

                if (node.isLeaf()) {
                    return _.reduce(this.options.leafLabelAttributes, function(m, attribute) {
                        return this.options.labelMakers[attribute].call(this, node, m, y);
                    }, x + 2, this) - x + 2;
                } else {
                    return this.paper.text(x - 3, y - 9, node.getParameters().name).attr({
                        'text-anchor': 'end',
                        'font-size': '8px'
                    }).getBBox().width;
                }
            },
            renderNode: function (node) {
                // If we are a HIGHLIGHTed node, draw a circle on us, and a path from the root
                if (_.find(node.getParameters().properties, function (p) {
                    return p.ref.split(":")[0] === "HIGHLIGHT";
                })) {
                    var path = node.traverseUp(function(node, memo) {
                        var our_x = this._x(node);
                        var our_y = this._y(node);
                        var parent_x;

                        if (node.getParent()) {
                            parent_x = this._x(node.getParent());
                        } else {
                            parent_x = (our_x - (node.length() * this.horizontalScalingFactor()));
                        }   
        
                        memo.push([parent_x, our_y], [our_x, our_y]);
                        return memo;
                    }, [], this);
                    this.options.highlightPathCallback.call(this, this.paper.path("M" + path.join("L")).attr(this.options.highlightStrokeStyle));
                    this.options.highlightNodeCallback.call(this, this.paper.circle(this._x(node), this._y(node), 3).attr(this.options.nodeStyle));
                }
            },
            // this maybe should go into the node class
            // it traverses up, and returns the most recent property
            // that the clade has
            mostRecentCladeProperty: function(node) {
                var clade_property = _.find(node.getParameters().properties, function (p) {
                    return p.applies_to === "clade";
                });

                if (clade_property) {
                    return clade_property
                }
                else if (node.getParent()) {
                    return this.mostRecentCladeProperty(node.getParent());
                }
            },
            '_cladeColors': {},
            nextCladeColor: function () {
                var i = 0;
                return function () {
                    var colors = [
                        'rgb(255,0,0)',
                        'rgb(0,0,255)',
                        'rgb(255,255,0)' ];
                    return colors[i++ % colors.length];
                }
            }(),
            edgeColor: function (node) {
                var property = this.mostRecentCladeProperty(node.getParent());
                if (property) {
                    if ( this._cladeColors[property.ref] === undefined ) {
                        this._cladeColors[property.ref] = this.nextCladeColor();
                    } 
                    return this._cladeColors[property.ref]
                }
                return 'rgb(0,0,0)';
            },
            renderEdgeToParent: function (node) {
                var our_x = this._x(node);
                var our_y = this._y(node);
                var path;
                if (node.getParent()) {
                    var parent_x = this._x(node.getParent());
                    var parent_y = this._y(node.getParent());
                    path = this.paper.path("M" + parent_x + "," + parent_y + "L" + parent_x + "," + our_y + "L" + our_x + "," + our_y).attr(this.options.edgeStyle).attr({stroke: this.options.edgeColor(node)});
                } else {
                    path = this.paper.path("M" + our_x + "," + our_y + "L" + (our_x - (node.length() * this.horizontalScalingFactor())) + "," + our_y).attr(this.options.edgeStyle).attr({stroke: this.options.edgeColor(node)});
                }
                this.options.pathCallback.call(this, path, node.length());
                return path;
            }
        },

        // Set up the widget
        _create: function () {
            // Should this not use [0]?
            this.paper = Raphael(this.element[0], $(this.element[0]).width(), $(this.element[0]).height());
            this._parse();
            this.refresh();
        },

        // Parses the PhyloXML from the server
        _parse: function () {
            function parse_clade(clade, parent) {
                // branch lengths of 0 do very bad things to drawing algorithms
                var branch_length = Math.max(clade["branch_length"] || 1e-10, 1e-10);
                var name = clade.name ? clade.name[0].Text : "";

                var my_data = {
                    parent: parent,
                    distance_to_parent: branch_length,
                    name: name
                };

                // I'm throwing all in with phyloXML here, and just using its schema.
                // For all of its faults, it makes like a lot easier at this moment
                if (clade.sequence) {
                    my_data["sequences"] = clade.sequence;
                }
                if (clade.taxonomy) {
                    my_data["taxonomies"] = clade.taxonomy;
                }
                if (clade.property) {
                    my_data["properties"] = clade.property;
                }

                var me = Node(my_data);
                _.each(clade.clade, function (child) {
                    parse_clade(child, me);
                });
                return me;
            }

            // Do we want to handle this messing up, or just vomit on the console?
            this.root = parse_clade($.xmlToJSON(this.options.data.xml).phylogeny[0].clade[0]);
        },

        // This is how much we have to scale the branch lengths in order to fit things inside our box
        horizontalScalingFactor: function (force) {
            // This is an expensive method, but we don't want to just memoize it, as
            // we can turn off and on tree label options
            if (this._cached_scaling_factor &&! force) {
                return this._cached_scaling_factor;
            }

            // We render the labels off the screen.
            // This isn't the best.  renderLabel should return
            // a set of the labels, and we calculate the bbox
            // from that.  This was quicker to code up, but will probably
            // change soon.

            this._cached_scaling_factor = 1; // We have to put in a fake value, as we have to render the labels.  This prevents infinite recursion

            var width = this.paper.width;

            // In general, we are just calculating the scaling factor that each leaf would need to
            // butt up against the right, and then storing the smallest
            var min_scaling_factor = _.reduce(this.root.leaves(), function (min, node) {
                var scaling_factor = ( width - this.options.renderLabel.call(this, node) - this.options.paddingLeft ) / node.distanceFromLeft();
                return Math.min(min, scaling_factor);
            }, Infinity, this);
            this.paper.clear();

            this._cached_scaling_factor = min_scaling_factor;
            return min_scaling_factor;
        },

        // "Node" doesn't know about the container, we do.  This takes in a node and returns coordinates for our container
        _x: function(node) {
            return node.distanceFromLeft() * this.horizontalScalingFactor() + this.options.paddingLeft;
        },
        _y: function(node) {
            // If we want to show the scale, we want to reserve a row on top for it
            var offset =  ( this.options.showScale ) ? 1.5 : 0.5;
            return this.options.rowHeight * (node.row() + offset);
        },
        refresh: function () {
            // If we aren't visible, sod off
            if ($(this.element[0]).is(":hidden")) {
                return
            }

            var width = $(this.element[0]).width();
            var height = ( this.options.showScale ) ? 
                ( this.root.leaves().length + 1 ) * ( this.options.rowHeight )  : // extra room if we need the scale
                this.root.leaves().length * ( this.options.rowHeight )  ;
            this.paper.setSize(width, height);

            // As you can imagine, this is the easy way.  We could scale our paths intead of redrawing them, which
            // would be much faster.  It would also be a bit of a pain to program, and we'd need to keep track of
            // the fake labels we make, and delete them one by one.
            this.paper.clear();

            // TODO detect width changes and only do this when it changes (or when the labels change)
            // This probably requires a richer event model
            // This is _just_ a performance increase
            this.horizontalScalingFactor(true);

            this.root.dft(function (node, p) {
                this.options.renderEdgeToParent.call(this, node);
                this.options.renderNode.call(this, node);
                this.options.renderLabel.call(this, node);
                return p;
            }, 0, this);

            if ( this.options.showScale ) {
                var  x1 = this.options.paddingLeft,
                     x2 = width / 4;

                var label = Math.round( ( x2 - x1 ) / this.horizontalScalingFactor() * 1000 ) / 1000;

                this.paper.path("M"+ x1 + ",0L" + x1 + "," + this.options.rowHeight).attr(this.options.edgeStyle);
                this.paper.path("M"+ x2 + ",0L" + x2 + "," + this.options.rowHeight).attr(this.options.edgeStyle);
                this.paper.path("M"+ x1 + "," + (this.options.rowHeight / 2) +"L" + x2 + "," + ( this.options.rowHeight / 2) ).attr(this.options.edgeStyle);
                this.paper.text(x1 + (x2 - x1) / 2, (this.options.rowHeight / 4), label);

            }
        },
        _setOption: function(option, value) {
            switch (option) {
                case "leafLabelAttributes":
                    this.options.leafLabelAttributes = value;
                    this.horizontalScalingFactor(true); // This line can be removed. It might be forcing recalculation, as the labels changed, even when not needed.
                    this.refresh();
                    break;
            }
        }
    });

    //
    // This is the Node, the basis of our tree.  This object tries to only know about
    // node things, and not worry about what happens to node.  It should be able to quickly
    // answer your questions about nodes, and relationships of nodes.  If you want to know
    // about how nodes act on the paper, than ask the container class. (e.g.: the location of
    // a node in pixels)
    var Node = function (params) {
            var children = []

            var that = {
                maxiumDistanceDown: function () {
                    return params.distance_to_parent + _.reduce(children, function (memo, child) {
                        return Math.max(memo, child.maxiumDistanceDown())
                    }, 0)
                },
                // Depth First Traversal
                // Many of the operations on the entire tree, or a clade, can be mapped onto a DFT.
                dft: function (lambda, memo, context) {
                    return lambda.call(context, this, _.reduce(children, function (memo, node) {
                        return node.dft(lambda, memo, context);
                    }, memo));
                },
                // Traverse Up
                // Apply the operation from yourself, up to the root
                traverseUp: function( lambda, memo, context) {
                    if ( this.getParent() ) {
                        return lambda.call(context, this, this.getParent().traverseUp(lambda, memo, context) );
                    } else {
                        return lambda.call(context, this, memo);
                    }
                },
                isLeaf: function () {
                    return (children.length == 0);
                },
                getChildren: function () {
                    return children
                },
                getParent: function () {
                    return params.parent;
                },
                getParameters: function () {
                    return params;
                },
                length: function () {
                    return params.distance_to_parent;
                },
                addChild: function (child) {
                    return children.push(child);
                },
                leaves: _.memoize(function () {
                    return this.dft(function (node, memo) {
                        if (node.isLeaf()) {
                            memo.push(node);
                        }
                        return memo;
                    }, [])
                }),
                root: _.memoize(function () {
                    if (params.parent) {
                        return params.parent.root();
                    }
                    return this;
                }),
                // This is the vertical offset value
                row: _.memoize(function () {
                    // If we are a leaf, we are an even number row, based on our order in the tree
                    if (this.isLeaf()) {
                        return _.indexOf(this.root().leaves(), this);
                    }
                    // If we are an internal node, we are the average of our childrens' rows
                    return _.reduce(children, function (m, c) {
                        return m + c.row()
                    }, 0) / children.length;
                }),
                distanceFromLeft: _.memoize(function () {
                    if (params.parent) {
                        return params.distance_to_parent + params.parent.distanceFromLeft();
                    }
                    return params.distance_to_parent;
                })
            };

            if (params.parent) {
                params.parent.addChild(that);
            }

            return that
        };

}(jQuery));



(function ($) {
    $.widget("ui.summarylist", {
        // These options will be used as defaults
        options: {
            visible: 3, //The number of items initially visible
            showMore: "Show More...",
            showFewer: "Show Fewer..."
        },
        _create: function() {
            if ( $(this.element).find('li').length > this.options.visible ) {
                $(this.element).find('li:gt('+ (this.options.visible - 1) +')').addClass("ui-summarylist-hidden-items").hide();
                $(this.element).append("<li><a href='javascript:;' class='ui-summarylist-show-more'>" + this.options.showMore + "</a></li>");
                $(this.element).append("<li style='display:none'><a href='javascript:;' class='ui-summarylist-show-fewer'>" + this.options.showFewer + "</a></li>");

                $(this.element).find(".ui-summarylist-show-more").click(function(ev) {
                    $(this).parent().siblings(':hidden').slideDown().end().slideUp();
                });
                $(this.element).find(".ui-summarylist-show-fewer").click(function(ev) {
                    $(this).parent().siblings('.ui-summarylist-hidden-items').slideUp().end().slideUp();
                    $(this).parent().parent().find('li .ui-summarylist-show-more').parent().slideDown();
                });
            }
        }
    });
}(jQuery));
