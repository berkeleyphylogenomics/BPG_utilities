function update_status() {
    $j.get(getString, {}, function(data) {
        $j('.bar').css('width', (data.status_id * (100/NUM_FATCAT_STAGES)).toString() + '%');
        $j('#status').text(data['status']);
        $j('#job-run-time').text(data['job-run-time']);
        if (data['stage1-families']) {
            $j(' #stage1-families ' ).text(data['stage1-families']);
            $j(' #stage1-families-table ').show();
        }
        if (data['stage2-families']) {
            $j(' #stage2-families ').text(data['stage2-families']);
            $j(' #stage2-families-table ').show();
        }
        if (data['done']) {
            $j('#status').text('Please wait while your results load.');
            location.reload(true);
        }
        else {setTimeout(update_status, 10000);} 
    });
}
   
$j(document).ready(function() {
    update_status();
});
