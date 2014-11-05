// File containing javascript for the intrepid webserver
// progress page.

function updateStatus() {
    $j.get(getString, {}, function(data) {
        // update the status and progress bar
        $j('.status-text').html(data['status']);
        $j('.bar').css('width', Math.floor(data['stage']*100/data['total_stages']).toString() + '%');
        // full length console
        if (data['full_length']) {
            // handle the full length update
            if ($j(' #full-length-console ').length == 0) {
                // we don't have the console yet, add it to the document
                $j(' #verbose-consoles ').append('<div class="console" id="full-length-console"><h5>Full Length Protein</h5>' +
                                    '<textarea rows="5" id="full-length-console-text" class="console-text"></textarea></div>');
                $j(' #verbose-consoles p ').remove();
            }
            $j(' #full-length-console-text ').val(data['full_length']);
        }        
        // process the domains
        for (var key in data) {
            if (key[0] == 'D') {
                // process the domain
                var domid = key.split("_")[1];
                if ($j(' #' + domid + '-console').length == 0) {
                    // we don't have the console yet, add it to the document
                    $j(' #verbose-consoles ').append('<div class="console" id="' + domid + '-console"><h5>' + 
                        data[key]['title'] + '</h5><textarea rows="5" id="' + domid + 
                        '-console-text" class="console-text"></textarea></div>');
                    $j(' #verbose-consoles p ').remove();
                }
                $j(' #' + domid + '-console-text').val(data[key]['status']);
            }
        }
        if (data['stage'] == data['total_stages']) {
            window.location = '/intrepid/' + data.id + '/';
        }
        else {setTimeout(updateStatus, progressAPIUpdateInterval);}
    });
}

$j(document).ready(function() {
    updateStatus();
    // bind click handlers for links that uncover divs
    $j(' .shjs ').click(function() {
        $j( $j(this).data().uncover ).toggle();
    });
});
