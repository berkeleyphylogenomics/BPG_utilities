// File containing javascript for the protein analysis webserver
// progress page.

function updateStatus() {
    $j.get(getString, {}, function(data) {
        if (data['domains_scanned']) {
            window.location = '/SAS/' + data.id + '/';
        }
        else {setTimeout(updateStatus, progressAPIUpdateInterval);}
    });
}

$j(document).ready(function() {
    updateStatus();
});
