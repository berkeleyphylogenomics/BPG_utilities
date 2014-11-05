function qtipStyle(obj) {
    return _.defaults(obj, {
        position: {
            my: 'top left',
            at: 'top right',
            target: 'mouse',
            adjust: {
                x: 15,
                y: 15
            },
            viewport: $j(window)
        },
        style: {
            classes: 'ui-tooltip-shadow ui-tooltip-bootstrap',
            tip: true
        }            
    });
}

// color functions
function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null
}

// jmol global variables
var jmolViewer;
var jmolInitializer;
var jmolSpin = false;

// jmol functions
function jmolReady(applet) {
}

function initializeJmolView() {
    var thisInitializationScript = 'select *;spacefill only; color structure; display :' + chainID + '; center :' + chainID + '; zoom 0; hide waters;' + getHighlightingScript();
    Jmol.script(jmolViewer, thisInitializationScript);
}

function getHighlightingScript() {
    var this_data_object;
    // unhighlight the residues
    Jmol.script(jmolViewer, 'select *; color ' + unrankedResidueColor);
    if ($j(' #structure-score-type-selector option:selected ').val() == 'discern') {
        this_data_object = discern_scores;
    }
    else if ($j(' #structure-score-type-selector option:selected ').val() == 'intrepid') {
        this_data_object = intrepid_scores;
    }

    // get the number of residues to highlight
    var loopCount = 0;
   
    var thisMax = this_data_object[0][1];
    var thisMin = 100;
    var thisLambda = parseFloat($j(' #lambda-value ').val());
    var thisScript = '';

    // find the minimum score that is actually reliable
    for (var i=0; i < this_data_object.length; i++) {
        // scores less than -100 are errors
        if ((this_data_object[i][1] > -100) && (this_data_object[i][1] < thisMin)) {thisMin=this_data_object[i][1];}
    }

    // calculate the x shift in terms of discern score
    var thisX0 = thisMax - ((thisMax - thisMin) * parseFloat($j(' #logisticx0 ').val()) / 100);

    var thisStartColors = hexToRgb($j(' #cold-color ').val());
    var thisEndColors = hexToRgb($j(' #hot-color ').val());

    // actually loop through them and highlight the residues
    for (var i=0; i < this_data_object.length; i++) {
        if (this_data_object[i][1] >= thisMin) {
            var thisScore = parseFloat(this_data_object[i][1]);
            // 
            var thisParam =  (1 / ( 1 + Math.exp( -1.0 * thisLambda * (thisScore - thisX0))));
            var thisRed = parseInt( Math.round( thisStartColors.r + thisParam * (thisEndColors.r - thisStartColors.r)));
            var thisBlue = parseInt( Math.round( thisStartColors.b + thisParam * (thisEndColors.b - thisStartColors.b)));
            var thisGreen = parseInt( Math.round( thisStartColors.g + thisParam * (thisEndColors.g - thisStartColors.g)));
            thisScript = thisScript + 'select ' + this_data_object[i][0].toString() + ';color [' +
                        thisRed.toString() + ',' + thisGreen.toString() + ',' + thisBlue.toString() + '];';
        }
    }
    return thisScript;
}

// binds the handlers between jquery interfaces and the jsmol/jmol applet
function bindHandlers() {
    // set the handler for the run button
    $j(' #execute-structure-command ').click(function() {
        if ($j(' #structure-command-line ').val().length != 0) {
            Jmol.script(jmolViewer, $j(' #structure-command-line ').val());
        } 
    });
    // do the toggle spin button
    $j(' #spin-structure ').click(function() {
        if (jmolSpin) {
            jmolSpin = false;
            Jmol.script(jmolViewer, 'spin off');
        }
        else {
            jmolSpin = true;
            Jmol.script(jmolViewer, 'spin on');    
        }
    });

    // do the display model stuff
    $j(' input:radio[name=display-model] ').click(function() {
        var dmv = $j(this).val();
        if (dmv == 'bs') {
            Jmol.script(jmolViewer, 'select *; cartoons off; spacefill 23%; wireframe 0.15');
        }
        else if (dmv == 'w') {
            Jmol.script(jmolViewer, 'select *; cartoons off; wireframe -0.1');
        }
        else if (dmv == 'c') {
            Jmol.script(jmolViewer, 'select protein or nucleic; cartoons only; set cartoonFancy true');
        }
        else if (dmv == 'sf') {
            Jmol.script(jmolViewer, 'select *; cartoons off; spacefill only');
        }        
    });

    // do the export
    $j(' #export-structure ').click(function() {
        Jmol.script(jmolViewer, 'write IMAGE PNGJ "jmol.png"');
    });

    // do the solvent accessible surface
    $j(' #toggle-solvent-surface ').click(function() {
        if (solventSurface === true) {
            Jmol.script(jmolViewer, 'isoSurface off');
            solventSurface = false;
        }
        else {
            Jmol.script(jmolViewer, 'isoSurface solvent');
            solventSurface = true;
        }
    });

    // handle the known functional sites
    $j(' #show-known-residues ').click(function() {
        Jmol.script(jmolViewer, 'select *; color ' + unrankedResidueColor);
        if (csaResidues.length>0) {Jmol.script(jmolViewer, 'select ' + csaResidues + ';color ' + knownResidueColor);}
    });

    // handle the highlight
    $j(' #highlight-top-k').click(function() {
        var thisScript = getHighlightingScript();
        Jmol.script(jmolViewer, thisScript);
    });
} 

// functions to load the hmtl5 based jsmol viewer
function loadJSmolViewer() {
    // hide the choose structure buttons
    $j(' #choose-structure ').hide();
    // show the other stuff
    $j(' #discern-structure-holder ').show();
    // start the jsmol initialization
    jmolInitializer = {
        addSelectionOptions: false,
        color: "#FFFFFF",
        debug: false,
        defaultModel: "",
        j2sPath: "/static/js/jsmol/j2s/",
        readyFunction: null,
        script: null,
        disableJ2SLoadMonitor: false,
        disableInitialConsole: false,
        readyFunction: jmolReady,
        allowjavascript: true,
        serverURL: "/static/js/jsmol/php/jsmol.php",
        src: null,
        use: "HTML5",
        width: "750",
        height: "750"
    };
    // don't put the jmol canvas on the page yet
    Jmol.setDocument(false);
    // actually make the applet
    Jmol.getApplet("jmolViewer", jmolInitializer);
    // put the object in the right place.
    if ($j(' #discern-structure-container ').length != 0) {
        $j(' #discern-structure-container ').html(Jmol.getAppletHtml(jmolViewer));
        // load the proper pdb structure
        Jmol.loadFile(jmolViewer, 'pdb/');
        initializeJmolView();
    }
    
    bindHandlers();
}

var pfam_map_args = {
    masterLength: masterSequenceLength,
    domainCallback: function(domain, dom_rect) {
        var content = domain.shortName + 
            " (" + domain.pfam_accession + ")<br>Alignment: " 
            + domain.alignment_from.toString() + " - " + 
            domain.alignment_to.toString() + "<br>E-value: " + 
            domain.evalue.toExponential() + "<br>" + domain.description;
        $j(dom_rect.node).qtip(qtipStyle({content: content}));
    },
    domainClick: function(domain) {
        // Here we should jump to the corresponding tab for this domain, it has been analyzed.
        // $j(' #domain-' + domain.id.toString()).click();
    }
};
pfam_map_args.domains = pfam_domains;

function startJalview(alignment,title,alwvar) {
      eval("var "+alwvar+" = document.JalviewLite.loadAlignment(alignment,title)");
}

function update_widgets() {
    if (!($j('#full-length').is(':hidden'))) {
        var w = $j('#full-length').width() - 147 - 10;
        // resize these widgets
        $j('#full-length table tr td:nth-child(2) div').width(w);
        $j('#pfam-match-widget').pfammap('refresh');
        $j('#query-line-widget').queryline('refresh');
    }
}

// debounce window resize
var debouncedResize = _.debounce(update_widgets, 500);

$j(document).ready(function() {
    $j( "#left-tabs" ).tabs({
        //'disabled': disabledArray,
        'create': update_widgets,
        'show':  update_widgets
    }).addClass('ui-tabs-vertical ui-helper-clearfix');
    $j( "#left-tabs li").removeClass('ui-corner-top').addClass('ui-corner-left');
    // Initialize widgets for the summary tab
    $j(' #pfam-match-widget ').pfammap(pfam_map_args);
    $j(' #query-line-widget ').queryline({masterLength:masterSequenceLength});
    // bind click handlers for links that click other links
    $j(' .jqls ').click(function() {
        $j( $j(this).data().linkto ).click();
    });

    // bind click handlers for links that uncover divs
    $j(' .shjs ').click(function() {
        $j( $j(this).data().uncover ).toggle();
    });

    // do the information tabs pane
    $j(' .intrepid-information-pane ').tabs({}).removeClass('ui-tabs-vertical ui-helper-clearfix');

    // Load JSMol view
    loadJSmolViewer();

    // debounce window resize
    window.addEventListener('resize', debouncedResize, false);

    // apply tooptips
    $j(".tip").qtip(qtipStyle({})); 
    
});
