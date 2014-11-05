// File containing javascript for the protein analysis webserver
// index page.

$j(document).ready(function() {
    // document ready function
    
    // Click handler for form submission
    $j(' #submit ').click(function() {
        // disable to submit button
        $j(this).attr('disabled', 'disabled');
        // Turn off previous alerts
        $j(' .protein-analysis-input-form tr td:nth-child(3) span ').html('');
        $j(' .protein-analysis-input-form td ').removeClass('pa-error-color');
        $j(' .protein-analysis-input-form input ').removeClass('form-error').removeClass('pa-error-color');
        $j(' .protein-analysis-input-form textarea ').removeClass('form-error').removeClass('pa-error-color');
        $j(' .protein-analysis-input-form .alert ').removeClass('alert-success').removeClass('alert-error').hide();
        // get the state of the checkboxes
        var domainsAnalyzed = '';
        $j(' .domains-scanned input[type=checkbox]:checked ').each( function(i) {domainsAnalyzed+=$j(this).val() + ' '});
        
        // post to the api
        $j.post('/api/SAS/', {
            'fasta': $j(' .fasta-input ').val(),
            'domains-to-scan': domainsAnalyzed,
            'minimum-domain-length': $j(' #minimum-domain-length-input ').val()
            }, function(data) {
                // callback function when api returns
                if (data.id) {
                    window.location = '/SAS/' + data.id + '/';
                }
                else {
                    if (data.status == 'error') {
                        if (data.type == '') {
                            // show the error
                        }
                    }
                    else {}
                }
            }
        );
            
        $j(this).removeAttr('disabled','disabled');
    }); 

    // Click handler for show/hide stuff on this page
    $j(' .jqsh ').click(function() {
        $j($j(this).data().toggle).toggle();
    });
    
});
