// Add the following to your HTML:
// <script type="text/javascript" src="http://path/to/archaeopteryx.js"></script>
//
// In method "openWin( dataurl, configfile ), replace "http://path/to/archaeopteryx_applets_jar_directory" with path of your archaeopteryx_applets.jar.
//
// Call method "openWin( dataurl, configfile )" with something like:
// <a href='#' onclick="openWin( 'http://path/to/treefile', 'http://path/to/config' )">launch Archaeopteryx</a>

function openArchaeopteryxWin( dataurl, configfile ) {
  
  aptx_window = open( "", "aptx_window", "width=300,height=150,status=no,toolbar=no,menubar=no,resizable=yes" );
  
  // open document for further output
  aptx_window.document.open();
  
  // create document
  aptx_window.document.writeln( "<HTML>" );
  aptx_window.document.writeln( "<HEAD>" );
  aptx_window.document.writeln( "<TITLE>Archaeopteryx Launchpad</TITLE>" );
  aptx_window.document.writeln( "</HEAD>" );
  aptx_window.document.writeln( "<BODY TEXT=\"#FFFFFF\" BGCOLOR=\"#000000\">" );
  aptx_window.document.writeln( "<FONT FACE=\"HELVETICA,ARIAL\">" );
  aptx_window.document.writeln( "<CENTER>" );
  aptx_window.document.writeln( "<B>Please do not close this window as long as you want to use Archaeopteryx.</B>" );
  aptx_window.document.write( "<APPLET ARCHIVE=\"archaeopteryx_applets.jar\"" );
  aptx_window.document.write( " CODEBASE=\"/static/java/archaeopteryx/\"" );
  aptx_window.document.write( " CODE=\"org.forester.archaeopteryx.ArchaeopteryxA.class\"" );
  aptx_window.document.write( " NAME=\"ArchaeopteryxA\"" );
  aptx_window.document.write( " WIDTH=\"220\" HEIGHT=\"60\"" );
  aptx_window.document.writeln( " ALT=\"ArchaeopteryxA is not working on your system (requires at least Sun Java 1.5)!\">" );
  aptx_window.document.writeln( "<PARAM NAME=\"url_of_tree_to_load\" VALUE=\"" + dataurl + "\">" );
  aptx_window.document.writeln( "<PARAM NAME=\"config_file\" VALUE=\"" + configfile + "\">" );
  aptx_window.document.writeln( "Your browser is completely ignoring the &lt;APPLET&gt; tag!" );
  aptx_window.document.writeln( "</APPLET>" );
  aptx_window.document.writeln( "</CENTER>" );
  aptx_window.document.writeln( "</BODY>" );
  aptx_window.document.writeln( "</HTML>" ); 
  
  // close the document - (not the window!)
  aptx_window.document.close();
}
  
  
