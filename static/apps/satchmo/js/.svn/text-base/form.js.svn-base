$(document).ready(function() {
    $('#submission_form').submit(function (e) {
        var msg = "";
        var blank_re = /^\s*$/;
        var textblank = blank_re.test(this.satchmo_fasta_fasta_text.value);
        var fileblank = blank_re.test(this.satchmo_fasta_fasta_file.value);
        var headless = this.satchmo_fasta_fasta_text.value.replace(/[\n\r]+/g,"\n").replace(/>[^\n]+\n/g,"").replace(/\s/g,"").toUpperCase();
        var nts = headless.replace(/[^ATCG]/g,"");

        if (!textblank && !fileblank) {
            //msg += "Please submit EITHER text input OR a file, not both.\n";
            msg += "Please submit text input.\n";
            this.satchmo_fasta_fasta_text.focus();
        } else if (textblank && fileblank) {
            msg += "Please submit either text input or a file.\n";
        } else if (!textblank && parseFloat(nts.length)/headless.length > 0.9) {
            msg += "Submission appears to be nucleotide sequences; only protein sequences are accepted.";
        }

        if (blank_re.test(this.email_to.value)) {
            msg += "Please provide an email address (for notification of job completion).\n";
        }
        /*
        if (blank_re.test(this.email_subject.value)) {
            this.email_subject.value = "BPG SATCHMO Results";
        }
        if (!blank_re.test(msg)) {
            alert(msg);
            return false;
        }
        */
        var bad_header = /^>.*?[\(\)]/m;
        if (bad_header.test(this.satchmo_fasta_fasta_text.value)) {
            this.satchmo_fasta_fasta_text.value = this.satchmo_fasta_fasta_text.value.replace(/\(/g,'%28').replace(/\)/g, '%29').replace(/:/g,'%3A').replace(/;/g,'%3B').replace(/\[/g,'%5B').replace(/\]/g,'%5D').replace(/,/g,'%2C');
        }
    });
    $('#example0').click(function() {
        f = $(this).closest('form')[0];
        $.get('/static/apps/satchmo/download/example0.fasta', function(data) {f.satchmo_fasta_fasta_text.value = data; });
        return false;
    });
    $('#example1').click(function() {
        f = $(this).closest('form')[0];
        $.get('/static/apps/satchmo/download/example1.fasta', function(data) {f.satchmo_fasta_fasta_text.value = data; });
        return false;
    });
    $('#formreset').click(function () {
        f = $(this).closest('form');
        f[0].reset();
        $(f).find('textarea').val('');
        return false;
    });
});

