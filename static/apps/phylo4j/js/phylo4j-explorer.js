/* Contains the javascript functions for the phylo4j explorer */
// variables for edge sets
var allPPI = ['AGREEMENT_SUBTREE_PPI','dip_ppi','PROTEIN_PROTEIN_INTERACTION'];

var number_loads_done = 0;
var thisNodePropertiesTable = [];
var thisRelationshipsTable = [];
var thisSelectionTable = [];
var thisGraphData = {};
var theseRelationships = [];
var nodeID = 0;

function getRelationshipTypes() {
    // this funtion will get the relationship types from the modal section
    // and return them
    var retVal = [];
    $j(' .edge-checkbox-container input ').each(function() {
        if($j(this).is(":checked")) {retVal.push($j(this).val());}
    });
    return retVal;
}

function displayPhylo4jExplorer() {
    // display the phylo4j explorer and turn off all the stuff
    $j(' .explorer-search-loading ').addClass('hidden');
    $j(' .explorer-loading-status ').addClass('hidden');
    $j(' .instruction-text ').removeClass('hidden');
    $j(' .explorer-display ').removeClass('hidden');
    $j(' .explorer-properties-main-table ').removeClass('hidden');
    $j(' #explorer-search-input ').prop('disabled', false);
    $j(' #explorer-search-button ').prop('disabled',false);
}

function loadPhylo4jExplorer(node) {
    // function gets all of the necessary data from the API to populate this page.
    // Get the properties for this node from the server
    $j.ajax('/api/phylo4j/explorer/' + node.toString() + '/properties/', {
        data: {},
        success: function(html) {
            $j(' .this-node-table ').html(html).removeClass('hidden');
            thisNodePropertiesTable = $j(' .this-node-table table ').dataTable({
                'dom':'lpt',
                'lengthMenu': [ [5,10,25,-1], [5,10,25,"All"] ],
                'pageLength':-1
            });
            number_loads_done++;
            if (number_loads_done == 3) {displayPhylo4jExplorer();}
        }
    });

    // Set up the datatable for the relationships
    $j.ajax('/api/phylo4j/explorer/' + node.toString() + '/relationships/', {
        data: {},
        success: function(data) {
            theseRelationships = data['relationships'];
            // disable the relevant checkboxes
            $j(' .edge-checkbox-container input ').each(function() {
                if(theseRelationships.indexOf($j(this).val())>=0) {$j(this).prop('disabled',false).prop('checked',true).parent().removeClass('hidden');}
                else {$j(this).prop('disabled',true).prop('checked',false).parent().addClass('hidden');}
            });
            $j(' .these-relations-table ').html(data['html']).removeClass('hidden');
            thisRelationshipsTable = $j(' .these-relations-table table ').dataTable({
                'dom':'lpt',
                'lengthMenu': [ [5,10,25,-1], [5,10,25,"All"] ],
                'pageLength':-1
            });
            number_loads_done++;
            if (number_loads_done == 3) {displayPhylo4jExplorer();}
        }
    });

    // Get the graph topology data from the server
    $j.ajax('/api/phylo4j/explorer/' + node.toString() + '/graph_data/', {
        data: {
            'depth': $j(' #graph-depth-control ').val(),
            'relationship_types': getRelationshipTypes()
        },
        success: function(data) {
            thisGraphData = data;
            if ('error' in data) {$j(' .explorer-graph-display ').html(data['error']);}
            else {
                // load the graph system
                pSys = arbor.ParticleSystem(particleSystemProperties);
                //pSys.parameters({gravity:true});
                pSys.renderer = explorerRenderer('explorer-graph-display');
                pSys.graft(data); 
            }
            number_loads_done++;
            if (number_loads_done == 3) {displayPhylo4jExplorer();}
        }
    });

    // Sets the explorer title
    $j(' .graph-title ').text('Exploring ' + node.toString()).removeClass('hidden');
}

function displaySearchError(m) {
// this function displays the text in m in the search error box, it also hides the other 2 things in its place
    $j(' .instruction-text ').addClass('hidden');
    $j(' .explorer-search-loading ').addClass('hidden');
    $j(' .explorer-loading-status ').addClass('hidden');
    $j(' .search-error ').text(m).removeClass('hidden');
    $j(' #explorer-search-button ').prop('disabled', false);
    $j(' #explorer-search-input ').prop('disabled', false);
}

// On document load
$j(document).ready(function() {
    // make the display section vertically resizable
    //$j(' .explorer-display ').resizable({
    //    handles: 's',
    //    minHeight: 500,
    //    maxHeight:1500
    //});
   
    // do tooltips for now...
    $j(' .tip ').tooltip({
        //track: true,
        //position: {my: 'left', at:'bottom', collision: 'fit'}
    });

    // assign dialog box to the modal section
    $j(' #edge-modal ').dialog({
        autoOpen:false,
        modal:true,
        buttons: {
            "Select All": function() {$j(' .edge-checkbox-container input:not([disabled]) ').prop('checked',true);},
            "Select all PPI": function() {$j(' .edge-checkbox-container input:not([disabled]) ').each(
                function() {
                    if($j.inArray($j(this).val(), allPPI) < 0) {$j(this).prop('checked',false);}
                    else {$j(this).prop('checked',true);}
                });
            },
            "Select None": function() {$j(' .edge-checkbox-container input:not([disabled]) ').prop('checked',false);},
            Close: function() {$j(this).dialog("close");}
        }
    });

    // pop open the modal window on the button click
    $j(' #customize-edges ').button().click(function() {$j(' #edge-modal ').dialog("open");});

    // refresh the graph data when the refresh button is clicked
    $j(' #refresh-graph ').button({
    }).click(function () {updateExplorerGraphDisplay();});

    // Bind click handler to search button
    $j(' #explorer-search-button ').click(function() {
        if($j.trim($j(' #explorer-search-input ').val()).length > 0) {
            // first turn off the alerts and disable the buttons
            $j(' #explorer-search-input ').prop('disabled', true);
            $j(' #explorer-search-button ').prop('disabled', true);
            $j(' .search-error ').addClass('hidden');
            $j(' .instruction-text ').addClass('hidden');
            $j(' .explorer-loading-status ').text('Determining graph entry point . . .').removeClass('hidden');
            $j(' .explorer-search-loading ').removeClass('hidden');
            // if there is something here, send it to the server to see if we can
            // find the node in the graph
            var thisQuery = $j.trim($j(' #explorer-search-input ').val());
            $j.ajax('/api/phylo4j/explorer/search/', {
                data: {'query':thisQuery},
                dataType: 'json',
                //error: displaySearchError('API error'),
                success: function(data) {
                    if('node_id' in data) {
                        $j(' .explorer-loading-status ').text('Loading data from server . . .');
                        nodeID = parseInt(data['node_id']);
                        loadPhylo4jExplorer(parseInt(data['node_id']));
                    }
                    else {
                        if('error' in data) {displaySearchError(data['message']);}
                        else {displaySearchError('API error');}
                    }
                }
            });
        }
    }); 
});
