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

// These are functions for the cluster display

var openRow = '';

function closeOpenClusterRow(oTable) {
    if ($j('.ortholog-cluster-row').length) {
        destroyClusterDetails(openRow, oTable);
        openRow = '';
    }
}    

function displayClusterDetails(data, nTr, oTable) {
    // displays the cluster details on the orhtologs tab
    var detailsRow = oTable.fnOpen(nTr, data['table'], 'ortholog-cluster-row');
    $j('.cluster-detail', detailsRow).slideDown();
    openRow = nTr;
    $j('.ortholog-cluster-row .tip').qtip(qtipStyle({}));
    return
}

function destroyClusterDetails(nTr, oTable) {
    // destroys the cluster details on the orthologs tab 
    $j('.ortholog-cluster-row .tip').qtip("destroy");
    $j('.cluster-detail', $j(nTr).next()[0]).slideUp( function () {
        oTable.fnClose(nTr);
    });
    return
}

function loadClusterRow(clusterID, nTr, oTable) {
    if (oTable.fnIsOpen(nTr) ) {
        // dispose of this row if its open already
        destroyClusterDetails(nTr, oTable);
    }
    else{
        // otherwise, load the cluster table from the server (html comes back though data)
        closeOpenClusterRow(oTable);
        $j.ajax({
            // url for data
            'url': '/api/fatcat2/' + jobID + '/cluster/' + clusterID + '/',
            'type': 'GET',
            'success': function(data) {displayClusterDetails(data, nTr, oTable);}
        });
    }
    return     
}

// end cluster display functions


var goDatatable = '';

var pfam_map_args = {
    masterLength: masterSequenceLength,
    domainCallback: function(domain, dom_rect) {
        var content = domain.shortName + 
            " (" + domain.pfam_accession + ")<br>Alignment: " 
            + domain.alignment_from.toString() + " - " + 
            domain.alignment_to.toString() + "<br>E-value: " + 
            domain.evalue.toExponential() + "<br>" + domain.description;
        $j(dom_rect.node).qtip(qtipStyle({content: content}));
    }
    //domainHeight: 140
};
pfam_map_args.domains = pfam_domains;

function startJalview(alignment,title,alwvar) {
      eval("var "+alwvar+" = document.JalviewLite.loadAlignment(alignment,title)");
}

function update_widgets() {
    if (!($j('#family-matches').is(':hidden'))) {
        // this is for firefox.
        var w = $j('#family-matches .tab-description').width() - 147 - 10;  // for padding and stuff
        // resize the family-matches widgets
        $j('#family-matches .hitmap-table tr td:nth-child(2) div').width(w);
        $j('#family-matches-pfam').pfammap("refresh");
        $j('#family-matches-query').queryline("refresh");
        $j('#family-matches-pfacts').alignmentmap("refresh");
    }
    else if (!($j('#enclosing-clades').is(':hidden'))) {
        var w = $j('#enclosing-clades .tab-description').width() - 147 - 10;  // for padding and stuff
        // resize the family-matches widgets
        $j('#enclosing-clades .hitmap-table tr td:nth-child(2) div').width(w);
        $j('#enclosing-clades-pfam').pfammap("refresh");
        $j('#enclosing-clades-query').queryline("refresh");
        $j('#enclosing-clades-pfacts').alignmentmap("refresh");
    }
    else if (!($j('#distant-clades').is(':hidden'))) {
        var w = $j('#distant-clades .tab-description').width() - 147 - 10;  // for padding and stuff
        // resize the family-matches widgets
        $j('#distant-clades .hitmap-table tr td:nth-child(2) div').width(w);
        $j('#distant-clades-pfam').pfammap("refresh");
        $j('#distant-clades-query').queryline("refresh");
        $j('#distant-clades-pfacts').alignmentmap("refresh");
    }
    else if (!($j('#results-summary').is(':hidden'))) {
        var w = $j('#results-summary .hitmap-table').width() - 147 - 10;
        $j('#results-summary .hitmap-table tr td:nth-child(2) div').width(w);
        $j('#results-summary-pfam').pfammap('refresh');
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
    $j('#results-summary-pfam').pfammap(pfam_map_args);
    // Initialize widgets for family matches tab
    $j('#family-matches-pfam').pfammap(pfam_map_args);
    $j('#family-matches-query').queryline({masterLength:masterSequenceLength});
    // Initialize pfam map for distant clades tab
    $j('#distant-clades-pfam').pfammap(pfam_map_args);
    $j('#distant-clades-query').queryline({masterLength:masterSequenceLength});

    // Initialize pfam map for enclosing clades tab
    $j('#enclosing-clades-pfam').pfammap(pfam_map_args);
    $j('#enclosing-clades-query').queryline({masterLength:masterSequenceLength});

    // Initialize summary list for go annotations
    $j('.go-summarylist').summarylist({});
    
    // Initialize datatable for family matches
    var familyMatchesTable = $j(' .family-matches-table ').dataTable( {
        "sDom": "t<'row'<'span3'i><'span5'p><'span2'r>>",
        "bServerSide": true,
        "aaSorting": [[ 4, 'asc' ]],
        "bDeferRender": true, // maybe better for ajax?
        "sAjaxSource": '/api/fatcat2/' + jobID + '/family_matches/',
        'bProcessing': true,
        "bFilter": false,  
        "sPaginationType": "bootstrap",
        "bAutoWidth": false,
        "oLanguage": {
            "sLengthMenu": "_MENU_ records per page"
        },
        "aoColumns": [
                {'sWidth': '100px', 'bSearchable': false},
                {'sWidth': '50px', 'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'sWidth': '60px', 'bSearchable': false},
                {'sWidth': '75px', 'bSearchable': false},
                {'sWidth': '60px', 'bSearchable': false},
                {'sWidth': '60px', 'bSearchable': false}
            ],
        'fnServerData': function (sSource, aoData, fnCallback) {
            $j.ajax( {
                'dataType': 'json',
                'type': 'GET',
                'url': sSource,
                'data': aoData,
                'success': function(result) {
                    // show processing on the heat map
                    $j('#family-matches-pfacts').alignmentmap("destroy");
                    $j('#family-matches .pfacts-loading').show();
                    $j('#family-matches-pfacts').alignmentmap({
                        'masterLength': masterSequenceLength,
                        'matches': result['heatmapData'],
                        'matchCallback': function(match, r) {
                            // tooltip matches
                            var content = "Family: " + match.family_name + "<br>Type: " + match.family_type +
                                    "<br>E-Value: " + match.evalue.toExponential();
                            $j(r.node).qtip(qtipStyle({content:content}));
                        },
                        'matchMouseover': function(match) { return function() { 
                            $j("#fam_" + match.id).addClass("highlighted-row");   
                        } },
                        'matchMouseout': function(match) { return function() { 
                            $j("#fam_" + match.id).removeClass("highlighted-row");   
                        } },
                        'matchMouseclick': function(match) { return function() {
                            window.open('/phylofacts/family/bpg0' + match.id + '/');
                        } }

                    });
                    $j('#family-matches .pfacts-loading').hide();
                    $j('#family-matches-pfacts').show();
                    $j('.family-matches-table tbody .tip ').qtip("destroy");
                    fnCallback(result);
                }
            });
        },
        'fnDrawCallback': function( oSettings ) {
            $j(' .family-matches-table tbody .tip ').qtip(qtipStyle({}));
            $j(' #family-matches-pfacts rect ').qtip('destroy');
        }
    });
    
    // Initialize datatable for enclosing clades
    var enclosingCladesTable = $j(' .enclosing-clades-table ').dataTable( {
        "sDom": "t<'row'<'span3'i><'span5'p><'span2'r>>",
        "bServerSide": true,
        "aaSorting": [[ 5, 'desc' ]],
        "bDeferRender": true, // maybe better for ajax?
        "sAjaxSource": '/api/fatcat2/' + jobID + '/enclosing_clades/',
        'bProcessing': true,
        "bFilter": false,  
        "sPaginationType": "bootstrap",
        "bAutoWidth": false,
        "oLanguage": {
            "sLengthMenu": "_MENU_ records per page"
        },
        "aoColumns": [
                {'sWidth': '80px', 'bSearchable': false},
                {'sWidth': '40px', 'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'sWidth': '60px', 'bSearchable': false},
                {'sWidth': '50px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '55px', 'bSearchable': false},
                {'bSearchable': false, 'bSortable':false, 'sWidth': '65px'}
            ],
        'fnServerData': function (sSource, aoData, fnCallback) {
            $j.ajax( {
                'dataType': 'json',
                'type': 'GET',
                'url': sSource,
                'data': aoData,
                'success': function(result) {
                    // show processing on the heat map
                    $j('#enclosing-clades-pfacts').alignmentmap("destroy");
                    $j('#enclosing-clades .pfacts-loading').show();
                    $j('#enclosing-clades-pfacts').alignmentmap({
                        'masterLength': masterSequenceLength,
                        'matches': result['heatmapData'],
                        'matchCallback': function(match, r) { //Add ToolTip
                            var content = "Family: " + match.family_name+ "<br>Subtree: " + match.node_name + "<br>E-Value: " + match.evalue.toExponential();
                            $j(r.node).qtip(qtipStyle({content: content}));
                        },
                        'matchMouseover': function(match) { return function() { 
                            $j("#ec_" + match.id).addClass("highlighted-row");   
                        } },
                        'matchMouseout': function(match) { return function() { 
                            $j("#ec_" + match.id).removeClass("highlighted-row");   
                        } },
                        'matchMouseclick': function(match) { return function() {
                            window.open('/phylofacts/tree_node_view/' + match.node_id.toString() + '/');
                        } }
                    });
                    $j('#enclosing-clades .pfacts-loading').hide();
                    $j('#enclosing-clades-pfacts').show();
                    $j('.enclosing-clades-table tbody .tip ').qtip("destroy");
                    fnCallback(result);
                }
            });
        },
        'fnDrawCallback': function( oSettings ) {
            $j(' .enclosing-clades-table tbody .tip ').qtip(qtipStyle({}));
            $j(' #enclosing-clades-pfacts rect ').qtip('destroy');
        }
    });

    // Initialize datatable for distant clades
    var distantCladesTable = $j(' .distant-clades-table ').dataTable( {
        "sDom": "t<'row'<'span3'i><'span5'p><'span2'r>>",
        "bServerSide": true,
        "aaSorting": [[ 6, 'desc' ]],
        "bDeferRender": true, // maybe better for ajax?
        "sAjaxSource": '/api/fatcat2/' + jobID + '/distant_clades/',
        'bProcessing': true,
        "bFilter": false,  
        "sPaginationType": "bootstrap",
        "bAutoWidth": false,
        "oLanguage": {
            "sLengthMenu": "_MENU_ records per page"
        },
        "aoColumns": [
                {'sWidth': '20px', 'bSearchable': false},
                {'sWidth': '80px', 'bSearchable': false},
                {'sWidth': '40px', 'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'sWidth': '60px', 'bSearchable': false},
                {'sWidth': '50px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '55px', 'bSearchable': false},
                {'bSearchable': false, 'bSortable':false, 'sWidth': '65px'}
            ],
        'fnServerData': function (sSource, aoData, fnCallback) {
            $j.ajax( {
                'dataType': 'json',
                'type': 'GET',
                'url': sSource,
                'data': aoData,
                'success': function(result) {
                    // show processing on the heat map
                    $j('#distant-clades-pfacts').alignmentmap("destroy");
                    $j('#distant-clades .pfacts-loading').show();
                    $j('#distant-clades-pfacts').alignmentmap({
                        'masterLength': masterSequenceLength,
                        'matches': result['heatmapData'],
                        'matchCallback': function(match, r) { //Add ToolTip
                            var content = "Family: " + match.family_name+ "<br>Subtree: " + match.node_name + "<br>E-Value: " + match.evalue.toExponential();
                            $j(r.node).qtip(qtipStyle({content: content}));
                        },
                        'matchMouseover': function(match) { return function() { 
                            $j("#dc_" + match.id).addClass("highlighted-row");   
                        } },
                        'matchMouseout': function(match) { return function() { 
                            $j("#dc_" + match.id).removeClass("highlighted-row");   
                        } },
                        'matchMouseclick': function(match) { return function() {
                            window.open('/phylofacts/tree_node_view/' + match.node_id.toString() + '/');
                        } }
                    });
                    $j('#distant-clades .pfacts-loading').hide();
                    $j('#distant-clades-pfacts').show();
                    $j('.distant-clades-table tbody .tip ').qtip("destroy");
                    fnCallback(result);
                }
            });
        },
        'fnDrawCallback': function( oSettings ) {
            $j(' .distant-clades-table tbody .tip ').qtip(qtipStyle({}));
            $j(' #distant-clades-pfacts rect ').qtip('destroy');
        }
    });
    
    // Initialize datatable for candidate orthologs
    var candidateOrthologsTable = $j(' .candidate-orthologs-table ').dataTable( {
        "sDom": "t<'row'<'span3'i><'span5'p><'span2'r>>",
        "bServerSide": true,
        "aaSorting": [[ 5, 'desc' ]],
        "bDeferRender": true, // maybe better for ajax?
        "sAjaxSource": '/api/fatcat2/' + jobID + '/candidate_orthologs/',
        'iDisplayLength': 50,
        'bProcessing': true,
        "bFilter": false,  
        "sPaginationType": "bootstrap",
        "bAutoWidth": false,
        "oLanguage": {
            "sLengthMenu": "_MENU_ records per page"
        },
        "aoColumns": [
                {'sWidth': '25px', 'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
            ],
        'fnServerData': function (sSource, aoData, fnCallback) {
            $j.ajax( {
                'dataType': 'json',
                'type': 'GET',
                'url': sSource,
                'data': aoData,
                'success': function(result) {
                    $j('.candidate-orthologs-table tbody .tip ').qtip("destroy");
                    fnCallback(result);
                }
            });
        },
        'fnDrawCallback': function( oSettings ) {
            $j(' .candidate-orthologs-table tbody .tip ').qtip(qtipStyle({}));
            // click handler for cluster display.
            $j(' .cluster-link ').click(function() {
                // this is really terrible, but couldn't get it to work any other way...
                loadClusterRow($j(this).attr('id'), $j(this)[0].parentNode.parentNode.parentNode.parentNode.parentNode.parentNode, candidateOrthologsTable);
            });
        }
    });
   
    // Initialize datatable for other sequence matches
    var otherSequenceMatchesTable = $j(' .other-sequence-matches-table ').dataTable( {
        "sDom": "t<'row'<'span3'i><'span5'p><'span2'r>>",
        "bServerSide": true,
        "aaSorting": [[ 5, 'desc' ]],
        "bDeferRender": true, // maybe better for ajax?
        "sAjaxSource": '/api/fatcat2/' + jobID + '/other_sequence_matches/',
        'iDisplayLength': 50,
        'bProcessing': true,
        "bFilter": false,  
        "sPaginationType": "bootstrap",
        "bAutoWidth": false,
        "oLanguage": {
            "sLengthMenu": "_MENU_ records per page"
        },
        "aoColumns": [
                {'sWidth': '25px', 'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
                {'sWidth': '45px', 'bSearchable': false},
            ],
        'fnServerData': function (sSource, aoData, fnCallback) {
            $j.ajax( {
                'dataType': 'json',
                'type': 'GET',
                'url': sSource,
                'data': aoData,
                'success': function(result) {
                    $j('.other-sequence-matches-table tbody .tip ').qtip("destroy");
                    fnCallback(result);
                }
            });
        },
        'fnDrawCallback': function( oSettings ) {
            $j(' .other-sequence-matches-table tbody .tip ').qtip(qtipStyle({}));
        }
    });
    // datatable for descriptions
    if ($j('#descriptions')) {
        var descriptionTable = $j(' #descriptions table ').dataTable( {
            'sDom': "t<'row'p>",
            'aaSorting': [[1, 'desc' ]],
            'iDisplayLength': 5,
            'bSort': false,
            'sPaginationType':'bootstrap',
            'bSearch': false,
            'bLengthChange': false,
            'bAutoWidth': false,
            "oLanguage": {
                "sLengthMenu": "_MENU_ records per page"
            },
            'aoColumns': [
                {'sWidth': '250px'},
                {'sWidth': '35px'}
            ]
        });
    }
 
    // datatable for genes
    if ($j('#genes')) {
        var descriptionTable = $j(' #genes table ').dataTable( {
            'sDom': "t<'row'p>",
            'aaSorting': [[1, 'desc' ]],
            'iDisplayLength': 5,
            'bSort': false,
            'bSearch': false,
            'sPaginationType': 'bootstrap',
            'bLengthChange': false,
            'bAutoWidth': false,
            "oLanguage": {
                "sLengthMenu": "_MENU_ records per page"
            },
            'aoColumns': [
                {'sWidth': '250px'},
                {'sWidth': '35px'}
            ]
        });
    }

    // click handlers
    // click handler for go annotation drill down on summary
    $j(' .go-annotation-view-hide-selector ').click(function() {
        if($j(this).hasClass('icon-zoom-in'))
        {
            var thisSelector = '#' + $j(this).data().accession.replace(':', '') + '-table';
            // is there a datatable open?  if so, close it...we don't want more than one open because
            // it seems better that way?
            if(goDatatable) {
                // get the table and destroy it.
                $j('#' + $j('.go-summarylist .icon-zoom-out').data().accession.replace(':','') + '-table table tbody .tip ').qtip("destroy");
                goDatatable.fnDestroy();
                $j('#' + $j('.go-summarylist .icon-zoom-out').data().accession.replace(':', '') + '-table').hide();
                $j('.go-summarylist .icon-zoom-out').removeClass('icon-zoom-out').addClass('icon-zoom-in').attr('title','View all orthologs with this GO annotation');
                goDatatable = '';
            }
            // show the table
            $j(thisSelector).show();
            goDatatable = $j(thisSelector + ' table').dataTable( {
                "sDom": "t<'row'<'span3'i><'span5'p><'span2'r>>",
                "bServerSide": true,
                "aaSorting": [[ 3, 'asc' ]],
                "bDeferRender": true, // maybe better for ajax?
                "sAjaxSource": '/api/fatcat2/' + jobID + '/annotations/go/' + $j(this).data().accession.split(':')[1] + '/',
                'bSort': true,
                'iDisplayLength': 5,
                'bProcessing': true,
                "bFilter": false,  
                "sPaginationType": "bootstrap",
                "bAutoWidth": false,
                "oLanguage": {
                    "sLengthMenu": "_MENU_ records per page"
                },
                "aoColumns": [
                        {'sWidth': '95px', 'bSearchable': false},
                        {'bSearchable': false},
                        {'bSearchable': false},
                        {'bSearchable': false}
                    ],
                'fnServerData': function (sSource, aoData, fnCallback) {
                    $j.ajax( {
                        'dataType': 'json',
                        'type': 'GET',
                        'url': sSource,
                        'data': aoData,
                        'success': function(result) {
                            $j(thisSelector + ' table tbody .tip ').qtip("destroy");
                            fnCallback(result);
                        }
                    });
                },
                'fnDrawCallback': function( oSettings ) {
                    $j(thisSelector + ' table tbody .tip ').qtip(qtipStyle({}));
                }
            });
            $j(this).removeClass('icon-zoom-in').addClass('icon-zoom-out').attr('title', 'Hide all orthologs with this GO annotation');
        }
        else {
            $j('#' + $j('.go-summarylist .icon-zoom-out').data().accession.replace(':','') + '-table table tbody .tip ').qtip("destroy");
            goDatatable.fnDestroy();
            $j('#' + $j('.go-summarylist .icon-zoom-out').data().accession.replace(':', '') + '-table').hide();
            $j('.go-summarylist .icon-zoom-out').removeClass('icon-zoom-out').addClass('icon-zoom-in').attr('title','View all orthologs with this GO annotation');
            goDatatable = '';
        }
    });

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
