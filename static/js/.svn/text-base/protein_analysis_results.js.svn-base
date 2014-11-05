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

function updateDomainStatus() {
    // This function is on a timer, it gets the status data from the api and 
    // updates the domain widgets accordingly
    $j.get('/api/SAS/' + jobID.toString() + '/', {'domains': 1}, function(data) {
        // callback function for status call
        for (var domain_id in data) {
            // is it a pfam domain ?
            if (data[domain_id]['type'] == 'pfam') {
                // update the pfam_domains variable for this domain
                for (var domain_dict in pfam_domains) {
                    if (domain_dict['id'] == domain_id) {
                        domain_dict['status_bit'] = data[domain_id]['status_bit'];
                        domain_dict['status'] = data[domain_id]['status'];
                        domain_dict['currentlyRunningJobStages'] = data[domain_id]['currentlyRunningJobStages'];
                        break;
                    }
                }   
            }
            else if (data[domain_id]['type'] == 'mda') {
                mda_domain['status_bit'] = data[domain_id]['status_bit'];
                mda_domain['status'] = data[domain_id]['status'];
                mda_domain['currentlyRunningJobStages'] = data[domain_id]['currentlyRunningJobStages'];
            }
        }
        update_widgets();
    });
    return {}
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
        if (domain.currentlyRunningJob == 'None') {
            
        {
        if (domain.status == 100) {
            // This domain has been analyzed,
        if (domain.status != 100) {
            openAnalyzeModal('Analyze ' + domain.shortName + '(' + domain.pfam_accession + ')', 
                domain.alignment_from.toString(), domain.alignment_to.toString(), domain.id,     
        }
        else {
            // Here we should jump to the corresponding tab for this domain, it has been analyzed.
            $j(' #domain-' + domain.id.toString()).click();
        }
    }
};
pfam_map_args.domains = pfam_domains;

var mda_map_args = {
    masterLength: masterSequenceLength,
    domainCallback: function(domain, dom_rect) {

    },
    domainClick: function(domain) {}
};
mda_map_args.domain = mda_domain;

function startJalview(alignment,title,alwvar) {
      eval("var "+alwvar+" = document.JalviewLite.loadAlignment(alignment,title)");
}

function update_widgets() {
    if (!($j('#protein-summary').is(':hidden'))) {
        var w = $j('#protein-summary').width() - 147 - 10;
        // resize these widgets
        $j('#protein-summary table tr td:nth-child(2) div').width(w);
        $j('#pfam-match-widget').pfammap('refresh');
        $j('#query-line-widget').queryline('refresh');
        $j('#mda-match-widget').mdamap('refresh');
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
    $j(' #mda-match-widget ').mdamap(mda_map_args); 
    // bind click handlers for links that click other links
    $j(' .jqls ').click(function() {
        $j( $j(this).data().linkto ).click();
    });

    // bind click handlers for links that uncover divs
    $j(' .shjs ').click(function() {
        $j( $j(this).data().uncover ).toggle();
    });

    // debounce window resize
    window.addEventListener('resize', debouncedResize, false);

    // apply tooptips
    $j(".tip").qtip(qtipStyle({})); 
    
});
