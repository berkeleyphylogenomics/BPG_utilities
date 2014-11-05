(function ($) {
    $.widget("ui.pairAliViewer", {
        options: {
            // text representation of the alignment
            alignment: "",
            // width of the widget
            paperWidth: 500,
            // alignment object
            alignmentObject: {},
            // identical alignment weight
            identicalWeight: 5,
            // non identical match weight
            matchWeight: 1
        },

        parseAlignment: function() {
            function getMatchPositions(seq1, seq2) {
                var matches = new Array();
                var matchString = "";
                for (i=0;i<seq1.length;i++) {
                    if ((seq1[i] != '-') && (seq2[i] != '-') && (seq1[i] != '.') && (seq2[i] != '.') && (seq1[i].toLowerCase() != seq1[i])) {
                        if (seq1[i] == seq2[i]) {matchString += seq1[i]; matches.push(that.options.identicalWeight);}
                        else {matchString += "+"; matches.push(that.options.matchWeight);}
                    }
                    else {matches.push(0); matchString+=" ";}
                }
                return {
                        'matches': matches,
                        'matchString': matchString
                       }
            }
            var that = this;
            var lines = this.options.alignment.split("\n");
            var headers = new Array();
            var sequences = new Array();
            var labels = new Array();
            var sequence = "";
            headers.push(lines[0]);
            for (var i=1;i<lines.length;i++) {   
                if ((lines[i].length) && (lines[i][0] == '>')) {
                    sequences.push(sequence);
                    sequence = "";
                    headers.push(lines[i]);
                }
                else {
                    sequence += lines[i];
                }
            }
            if (sequence) {sequences.push(sequence);}
            for (i=0;i<sequences.length;i++) {sequences[i] = sequences[i].replace(/[\s\n]/gm,'');}
            for (i=0;i<headers.length;i++) {labels.push(headers[i].replace(/>/g,'').split(" ")[0]);}
            var m = getMatchPositions(sequences[0], sequences[1]);
            var ret =  { 
                        'headers': headers,
                        'labels': labels,
                        'sequences': sequences,
                        'matchString': m.matchString,
                        'matches': m.matches
                       };
            console.log(ret);
            return ret 
        },

        _create: function () {
            //this.paper = Raphael(this.element[0], $(this.element[0]).width(), $(this.element[0]).height());
            //this.refresh();
            this.options.alignmentObject = this.parseAlignment();
        },

        refresh: function () {
            if ($(this.element[0]).is(":hidden")) {
                return
            }

            
        }
    });
}(jQuery));
