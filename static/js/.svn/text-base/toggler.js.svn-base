$(document).ready(function () {
    $('.toggle_text').show().children('.toggle_hidden').hide();
    $('.toggler').css('outline', 'none').each(function (i) {
        $('#'+this.id.split('_')[1]).addClass('togglee');
    }).click(function (e) {
        $(this).toggleClass('toggled').children('.toggle_text').children('span').toggle();
        $('#'+this.id.split('_')[1]).toggle();
    }).each(function (i) {
        if (! $(this).hasClass('toggle_defaultOn')) {
            $(this).click();
        }
    });
});
