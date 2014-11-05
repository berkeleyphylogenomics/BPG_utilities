      function annotate ( type, id ) {
         var annotateURL = "http://phylogenomics.berkeley.edu/annotate/form.php?type=" + type + "&id=" + id;
         annotateWindow = window.open( annotateURL, "AddAnnotation",
                                "scrollbars=yes, menubar=yes, resizable=yes, "
                                + "height=400, width=550" );
         annotateWindow.focus();
      }



      function show_comment ( type, id ) {
         var annotateURL = "http://phylogenomics.berkeley.edu/annotate/show_comment.php?type=" + type + "&id=" + id;
         annotateWindow = window.open( annotateURL, "AddAnnotation",
                                "scrollbars=yes, menubar=yes, resizable=yes, "
                                + "height=400, width=550" );
         annotateWindow.focus();
      }

