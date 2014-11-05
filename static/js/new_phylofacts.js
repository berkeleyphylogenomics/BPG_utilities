// Eventually, we should have a single global variable (say, BPG) and 
// all functions should be attributes of that global variable.

// This function makes several input boxes on index.html visible.

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
            rowHeight: 4,
            rowSpacing: 3,
            masterSequenceHeight: 12,
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
            this.paper.rect(0, 0, width, this.options.masterSequenceHeight);
            this.paper.text(10, Math.floor(this.options.masterSequenceHeight/2), "1");
            this.paper.text(width - 20, Math.floor(this.options.masterSequenceHeight/2), this.options.masterLength);

            _.each(this.options.matches, function (match, index) {
                var match_rect = this.paper.rect(
                match.alignment_match_from * scale, this.options.masterSequenceHeight + this.options.rowSpacing + (index) * (this.options.rowHeight + this.options.rowSpacing), (match.alignment_match_to - match.alignment_match_from) * scale, this.options.rowHeight).attr({
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
