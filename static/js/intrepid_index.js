// File containing javascript for the protein analysis webserver
// index page.

function postToAPI(fileStr) {
    $j.post('/api/intrepid/', {
        'email': $j(' .email-input ').val(),
        'minimum-column-coverage': $j(' #minimum-column-coverage ').val(),
        'minimum-pfam-domain-length': $j(' #minimum-pfam-domain-length ').val(),
        'homolog-iterations': $j(' #homolog-iterations ').val(),
        'pwid-query-homolog': $j(' #pwid-query-homolog ').val(),
        'sequence-database': (parseInt($j(' input[name="seqdb"]:checked ').val()) + 
                             (parseInt($j(' input[name="blastv"]:checked ').val()) *
                              $j(' input[name="seqdb"] ').length)),
        'homolog-coverage': $j(' #minimum-homolog-coverage ').val(),
        'maximum-homologs': $j(' #maximum-homologs ').val(),
        'save-program-outputs': $j(' #save-program-outputs ').is(':checked'),
        'pdb': $j(' .pdb-input ').val().replace(/ /g,''),
        'pdb-file': fileStr,    
    }, function(data) {
            // callback function when api returns
            if (data.id) {
                window.location = '/intrepid/' + data.id + '/';
            }
            else {
                if (data.status == 'error') {
                    if (data.type == 'fasta') {
                    }
                    else if (data.type == 'email') {
                        // show the error
                        $j(' .email-error ').addClass('form-error').addClass('intrepid-error-color').text(data.message).show();
                        $j(' .email-input ').addClass('form-error');
                    }
                    else if (data.type == 'pdb') {
                        $j(' .pdb-error ').addClass('intrepid-error-color').text(data.message).show();
                        $j(' .pdb-input ').addClass('form-error');   
                    }
                    else {
                        $j(' #job-create-message ').addClass('alert-error').text(data.message).show();
                    }
                }
                else {}
            }
        }
    );
}


$j(document).ready(function() {
    // document ready function
    
    // Click handler for form submission
    $j(' #submit ').click(function() {
        // disable to submit button
        $j(this).attr('disabled', 'disabled');
        // Turn off previous alerts
        $j(' .intrepid-input-form tr td:nth-child(3) span ').html('');
        $j(' .intrepid-input-form td ').removeClass('intrepid-error-color');
        $j(' .pdb-error ').html('').removeClass('intrepid-error-color').hide();
        $j(' .intrepid-input-form input ').removeClass('form-error').removeClass('intrepid-error-color');
        $j(' .intrepid-input-form textarea ').removeClass('form-error').removeClass('intrepid-error-color');
        $j(' .intrepid-input-form .alert ').removeClass('alert-success').removeClass('alert-error').hide();
        
        // get the file stuff
        var file_string = ''
        if ($j(' .pdb-file-input ')[0].files.length > 0) {
            var file = $j(' .pdb-file-input ')[0].files[0];
            var reader = new FileReader();
            reader.readAsText(file);
            reader.onloadstart = function() {
                $j(' .file-loading-text ').html('Starting to upload file ' + file);
            };

            reader.onloadend = function() {
                $j(' .file-loading-text ').html('');
            };

            reader.onload = function(e) {
                file_string = e.target.result
                postToAPI(file_string);
            };

            reader.onerror = function() {
                $j(' .pdb-error ').addClass('intrepid-error-color').text('Error reading file ' + file).show();
                $j(' .pdb-input ').addClass('form-error');
            };

        }
        else {
            postToAPI(file_string);
        }
        
        $j(this).removeAttr('disabled','disabled');
    }); 

    $j(' #auto-populate-pdb ').click(function() {
        $j(' .pdb-input ').val('12as:A');
    });

    // Click handler for show/hide stuff on this page
    $j(' .shjs ').click(function() {
        $j($j(this).data().uncover).toggle();
    });
    
});
