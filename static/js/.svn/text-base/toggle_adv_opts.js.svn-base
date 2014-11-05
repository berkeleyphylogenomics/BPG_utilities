/*
 * Toggle content as an element is clicked.
 */
var contentID   = 'adv_opts_table';
var buttonID    = 'adv_opts_button';
function toggleOptions() {
    var button      = document.getElementById(buttonID)
    var buttonText  = button.innerHTML;
    if (buttonText.match(/: Show$/)) {
        document.getElementById(contentID).style.display = 'block';
        button.innerHTML          = buttonText.replace(/Show$/g,'Hide');
        button.style.borderBottom = 'none';
    } else if (buttonText.match(/: Hide$/)) {
        document.getElementById(contentID).style.display = 'none';
        button.innerHTML          = buttonText.replace(/Hide$/g,'Show');
        button.style.borderBottom = 'solid 2px';
    }
}

window.onload=function() {

    // Put a new anchor inside the 'button'-labeled element
    var outerButton = document.getElementById(buttonID);
    outerButton.removeAttribute('id');
    var buttonText = outerButton.innerHTML + ': Show';
    outerButton.innerHTML = '';
    var newAnchor = document.createElement('a');
    newAnchor.setAttribute('href','#');
    newAnchor.setAttribute('id',buttonID);
    newAnchor.innerHTML             = buttonText;
    newAnchor.style.borderBottom    = 'solid 2px';
    newAnchor.style.width           = buttonText.length+'em';
    newAnchor.addEventListener('click',toggleOptions,false);
    outerButton.appendChild(newAnchor);

    // Style the content appropriately and hide
    var contentStyle = document.getElementById('adv_opts_table').style;
    contentStyle.display            = 'none';
    
    // Styling that should probably be moved to a CSS file
    newAnchor.style.display         = 'block';
    newAnchor.style.backgroundColor = '#DDDDDD';
    newAnchor.style.borderTop       = 'solid 2px';
    newAnchor.style.textDecoration  = 'none';
    newAnchor.style.color           = 'black';
    newAnchor.style.position        = 'relative';
    newAnchor.style.top             = '6px';
    newAnchor.hideFocus             = true;
    newAnchor.style.outline         = 'none';
    contentStyle.borderTop          = 'solid 2px';
    contentStyle.borderBottom       = 'solid 2px';
    contentStyle.backgroundColor    = '#DDDDDD';

    // Start with advanced options shown if there was an error
    var tds = document.getElementById(contentID).getElementsByTagName('td');

    for (var key in tds) {
        if (tds[key].className == 'adv_opts_err') {
            toggleOptions();
            break;
        }
    }
}
