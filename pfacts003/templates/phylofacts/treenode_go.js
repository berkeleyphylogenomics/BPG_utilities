<script>
$j(document).ready(function() {
        $j(".go-sequence-table").each(function(idx, el) {
            $j(el).dataTable( {
            "sDom": "<'row'<'span4'l><'span4'f>r>t<'row'<'span4'i><'span6'p>>",
            "sPaginationType": "bootstrap",
            "oLanguage": {
            "sLengthMenu": "_MENU_ records per page"
            }
            } );
        });
        $j(".go-table-toggle").click(function () {
                $j(this).parent().find(".go-sequence-table-container").toggle();
            });
        } );
</script>

