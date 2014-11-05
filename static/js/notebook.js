/*
 * Toggle content as an element is clicked.
 */

$(document).ready(function(){
    $('.notebook').each(function (notebook_i) {
        // Add a notebookTabBar
        var notebookTabBar = $(document.createElement('div')).addClass('notebookTabBar').prependTo(this)[0];
        // Move the notebookLabel into the notebookTabBar
        $(this).children('.notebookLabel').appendTo(notebookTabBar);
        // Hide notebookPages, move notebookTabs to the notebookTabBar,
        //  set initial conditions and events
        $(this).children('.notebookPage').toggle().siblings('.notebookTab').addClass('notebookUnselected').appendTo(notebookTabBar).click(function (ev) {
            $(this).siblings('.notebookTab').addClass('notebookUnselected').each(function (i) {
                $('#'+this.id.split('_')[1]).hide();
            });
            if ($(this).parent().parent().hasClass('notebookAlwaysOn')) {
                alert('eh');
                $('#'+this.id.split('_')[1]).show();
                $(this).removeClass('notebookUnselected');
            } else {
                $('#'+this.id.split('_')[1]).toggle();
                $(this).toggleClass('notebookUnselected');
            }
            ev.preventDefault();
        });
        // Open the page that should be on by default
        $(notebookTabBar).children('.notebookDefaultOn').click();
    });
});

/*
$(document).ready(function(){
    // give each notebook an ID
    $('.notebook').attr('id', function(notebook_i) {
        return 'notebook_'+notebook_i
    // Prepare the notebook
    }).each(function(notebook_i) {
        // Add the notebook tabs div
        $(document.createElement('div')).addClass('notebook_tabs').prependTo(this);
        // Add the notebook label
        $(this).children('.notebook_title').appendTo($(this).children('.notebook_tabs'));
        // Convert page labels to anchor tags in the notebook tabs
        $(this).children('.notebook_page').children('.notebook_tab').each(function() {
            $(document.createElement('a')).attr('href','#').attr('hidefocus','true').css('outline','none').html(this.innerHTML).addClass('notebook_tab').appendTo($(this).parent().parent().children('.notebook_tabs'));
        }).remove();


        $(this).children('.notebook_page').attr('id', function(page_i){
            return 'notebook_'+notebook_i+'_page_'+page_i;
        }).hide().parent().children('.notebook_tabs').children('.notebook_tab').addClass('notebook_tab_off').attr('id', function(page_i){
            return 'notebook_'+notebook_i+'_tab_'+page_i;
        }).click(function (e) {
            var notebook_indices = $(this).attr('id').split('_');
            var page_target = $('#notebook_'+notebook_indices[1]+'_page_'+notebook_indices[3]);
            if (!$(this).hasClass('notebook_tab_off')) {
                if(!$(page_target).parent().hasClass('notebook_always_on')) {
                    $(this).addClass('notebook_tab_off');
                    $(page_target).hide();
                }
            } else {
                $(this).parent().children('.notebook_tab').addClass('notebook_tab_off');
                $(page_target).parent().children('.notebook_page').hide();
                $(this).removeClass('notebook_tab_off');
                $(page_target).show();
            }

            e.preventDefault();
        });
        if ($(this).hasClass('notebook_always_on')) {
            $(this).children('.notebook_tabs').children('.notebook_tab:first').removeClass('notebook_tab_off');
            $(this).children('.notebook_page:first').show();
        }
    });
});*/

/*
var area_id = 'view_divs';
function toggleDiv() {
    var buttonText = this.id.split('_',1)[0];
    var pane = document.getElementById(area_id);
    var buttons = document.getElementById('button_row').getElementsByTagName('a');
    for ( key in buttons ) {
        var button = buttons[key];
        if (! button.id.match(/^\w+?_button$/)) continue;
        var otherButtonText = button.id.split('_',1)[0];
        if (otherButtonText == buttonText) {
            document.getElementById(otherButtonText+'_block').style.display = 'block';
            button.style.bottom = "-2px";
        } else {
            document.getElementById(otherButtonText+'_block').style.display = 'none';
            button.style.bottom = "0";
        }
    }
}

window.onload=function() {
    var pane = document.getElementById(area_id);
    var divs = pane.firstChild.getElementsByTagName('div');
    var boxes = new Array();
    for ( key in divs ) {
        var box = divs[key];
        if ((box.id != undefined) && (box.id.match(/^\w+?_block$/))) boxes.push(box);
    }
    var buttonRow = document.createElement('div');
    var divsDivStyle = pane.getElementsByTagName('div')[0].style;
    buttonRow.id = 'button_row';
    var label = document.createElement('p');
    label.innerHTML = '&nbsp;ALIGNMENTS:&nbsp;';
    label.style.display = 'inline';
    buttonRow.appendChild(label);
    divsDivStyle.borderTop = 'solid 2px #888888';
    pane.insertBefore(buttonRow, pane.firstChild);

    for ( key in boxes ) {
        var box = boxes[key];
        if (box.id == undefined) continue;
        box.style.display = 'none';
        var button = document.createElement('a');
        var buttonText = box.id.split('_',1)[0];
        button.innerHTML = '&nbsp;&nbsp;'+buttonText+'&nbsp;&nbsp;';
        button.setAttribute('id', buttonText+'_button');
        button.setAttribute('href', '#');
        button.hideFocus = true;
        
        buttonStyle = button.style;
        buttonStyle.textDecoration = 'none';
        buttonStyle.fontWeight = 'bold';
        buttonStyle.color = '#333333';
        buttonStyle.textAlign = 'center';
        buttonStyle.backgroundColor = '#bbbbbb';
        buttonStyle.borderTop = 'solid 2px #888888';
        buttonStyle.position = "relative";
        buttonStyle.display = 'inline';
        buttonStyle.marginRight = '2px';
        buttonStyle.outline = 'none';

        button.addEventListener('click', toggleDiv, false);
        buttonRow.appendChild(button);
    }
    boxes[0].style.display = "block";
    buttonRow.getElementsByTagName('a')[0].style.bottom = "-2px";
}
*/
