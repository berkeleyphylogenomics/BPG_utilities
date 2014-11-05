$(document).ready(function() {
    $('.expandable').children('.expand').siblings().toggle().parent().children('.expand, .collapse').click(function () {
        $(this).parent().children().toggle();
        return false;
    });
});
