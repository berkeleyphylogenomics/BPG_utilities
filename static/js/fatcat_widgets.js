(function ($) {
    $.widget("ui.pfammap", {
        options: {
            // sequence length
            masterLength: 1,
            // domain color set
            colorSet: ['#08d9ec', '#8cf504', '#fe0000', '#06f213', '#0a4de8', '#097be9', 
                        '#09a9ea', '#08d9ec', '#08edd1', '#07eea2', '#07ef73', 
                        '#06f143', '#29f305', '#5af405', '#8cf504', '#bef703', 
                        '#f9cd02', '#fa9b02', '#fc6801', '#fd3401' ],
            domainHeight: 20, 
            domains: [ ],
            domainBorderRadius: 8,
            queryBaseColor: "rgb(208,208,208)",
            queryLineHeight: 4,
            fontSize: 10, 
            domainClick: function (match) {
                return function() {}
            },
            domainMouseover: function (match) {
                return function() {}
            },
            domainMouseup: function (match) {
                return function() {}
            },
            domainMouseout: function (match) {
                return function() {}
            },
            domainCallback: function(match, node) {}
        },

        _create: function () {
            this.paper = Raphael(this.element[0], $(this.element[0]).width(), $(this.element[0]).height());
            this.refresh();
        },

        refresh: function () {
            if ($(this.element[0]).is(":hidden")) {
                return
            }
            this.paper.clear();
            var width = $(this.element[0]).width() - 10; 
            var height = Math.max(this.options.domainHeight, this.options.queryLineHeight);
            var scale = width / this.options.masterLength; 
            var queryY = Math.floor(height/2) - Math.floor(this.options.queryLineHeight/2);
            var domainY = Math.floor(height/2) - Math.floor(this.options.domainHeight/2);

            this.paper.setSize(width + 10, height);
            // Draw the query guide line
            this.paper.rect(0, queryY, width, this.options.queryLineHeight).attr({ fill: this.options.queryBaseColor, stroke: 'none' });

            var domains_encountered = {}; 
            var colorIndex = 0;

            for(i=0; i<this.options.domains.length; i++) {
                // do the ith domain
                if (this.options.domains[i].pfam_accession in domains_encountered) {
                    currentColor = domains_encountered[this.options.domains[i].pfam_accession];
                }
                else {
                    currentColor = this.options.colorSet[colorIndex % this.options.colorSet.length];
                    domains_encountered[this.options.domains[i].pfam_accession] = currentColor;
                    colorIndex += 1;
                }
                // draw the rectangle of the appropriate color
                var rect = this.paper.rect(this.options.domains[i].alignment_from * scale,
                    domainY, (this.options.domains[i].alignment_to
                    - this.options.domains[i].alignment_from) * scale,
                    this.options.domainHeight, this.options.domainBorderRadius).attr({
                        fill: currentColor, 'stroke': 'none', 'cursor':'default'
                    }).toFront();
                var text = this.paper.text(-100,-100,this.options.domains[i].shortName).attr({
                    'font-size': this.options.fontSize, 'text-anchor': 'start', 'cursor': 'default'
                });
                if (((rect.getBBox().width - text.getBBox().width) > 4) &&
                    ((rect.getBBox().height - text.getBBox().height) > 4))
                {
                    text.attr({'x': Math.floor(rect.getBBox().x + Math.floor((rect.getBBox().width - text.getBBox().width)/2)),
                               'y': Math.floor(height/2)}).toFront();
                }
                else {
                    text.remove();
                }
                // this is the rectangle for tooltips
                var rect2 = this.paper.rect(this.options.domains[i].alignment_from * scale,
                    domainY, (this.options.domains[i].alignment_to
                    - this.options.domains[i].alignment_from) * scale,
                    this.options.domainHeight, this.options.domainBorderRadius).attr({
                        'fill-opacity': 0, 'stroke-opacity': 0, 'cursor':'default',
                        'fill': 'black', 'stroke':'black'
                    }).toFront();
                
                rect.mouseover(this.options.domainMouseover(this.options.domains[i]));
                rect.mouseup(this.options.domainMouseup(this.options.domains[i]));
                rect.mouseout(this.options.domainMouseout(this.options.domains[i]));
                rect.click(this.options.domainClick(this.options.domains[i]));
                this.options.domainCallback(this.options.domains[i], rect2);
            }
        }
    });
}(jQuery));

// Query Line
//
// Draws the little query line box with start and end positions
(function ($) {
    $.widget("ui.queryline", {
        // These options will be used as defaults
        options: {
            // Length of the master sequence
            masterLength: 1,
            // width of the master sequence
            queryLineHeight: 12,
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
            var height = this.options.queryLineHeight;
            var scale = width / this.options.masterLength;

            this.paper.setSize(width + 10, height+2);

            // Master Sequence
            this.paper.rect(0, 1, width, this.options.queryLineHeight);
            this.paper.text(10, Math.floor(this.options.queryLineHeight/2), "1");
            this.paper.text(width - 20, Math.floor(this.options.queryLineHeight/2), this.options.masterLength);

        }
    });
}(jQuery));

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
            var height = this.options.matches.length * (this.options.rowHeight + this.options.rowSpacing); // matches
            //this.options.rowHeight * 4; //header and bottom area
            var scale = width / this.options.masterLength;

            this.paper.setSize(width + 10, height);

            _.each(this.options.matches, function (match, index) {
                var match_rect = this.paper.rect(
                match.alignment_match_from * scale, (index) * (this.options.rowHeight + this.options.rowSpacing), (match.alignment_match_to - match.alignment_match_from) * scale, this.options.rowHeight).attr({
                    fill: "rgb(" + this.options.ramp(match.evalue).join(',') + ")",
                    'stroke': 'none', 'cursor':'pointer'
                });
                match_rect.mouseover(this.options.matchMouseover(match));
                match_rect.mouseup(this.options.matchMouseup(match));
                match_rect.mouseout(this.options.matchMouseout(match));
                match_rect.click(this.options.matchMouseclick(match));
                this.options.matchCallback(match, match_rect);
            }, this);
        },

        destroy: function() {
            this.paper.remove();
            $.Widget.prototype.destroy.call( this );
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
                    $(this).parent().siblings('li:hidden').slideDown().end().slideUp();
                });
                $(this.element).find(".ui-summarylist-show-fewer").click(function(ev) {
                    $(this).parent().siblings('.ui-summarylist-hidden-items').slideUp().end().slideUp();
                    $(this).parent().parent().find('li .ui-summarylist-show-more').parent().slideDown();
                });
            }
        }
    });
}(jQuery));
