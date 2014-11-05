Clazz.declarePackage ("JSV.export");
Clazz.load (["JSV.api.JSVExporter", "JSV.export.FormContext"], "JSV.export.FormExporter", ["java.io.IOException", "JSV.common.JSVFileManager", "J.util.Logger"], function () {
c$ = Clazz.decorateAsClass (function () {
this.context = null;
this.errMsg = null;
this.currentTime = null;
this.out = null;
this.viewer = null;
Clazz.instantialize (this, arguments);
}, JSV["export"], "FormExporter", null, JSV.api.JSVExporter);
Clazz.prepareFields (c$, function () {
this.context =  new JSV["export"].FormContext ();
});
$_M(c$, "initForm", 
function (viewer, out) {
this.viewer = viewer;
this.out = out;
this.currentTime = viewer.apiPlatform.getDateFormat (false);
}, "JSV.common.JSViewer,JU.OC");
$_M(c$, "writeForm", 
function (templateFile) {
var error =  new Array (1);
var template = JSV.common.JSVFileManager.getResourceString (this, "resources/" + templateFile, error);
if (template == null) {
J.util.Logger.error (error[0]);
return error[0];
}this.errMsg = this.context.setTemplate (template);
if (this.errMsg != null) {
J.util.Logger.error (this.errMsg);
return this.errMsg;
}this.errMsg = this.context.merge (this.out);
if (this.out == null) return this.errMsg;
if (this.errMsg != null) {
J.util.Logger.error (this.errMsg);
throw  new java.io.IOException (this.errMsg);
}this.out.closeChannel ();
return "OK " + this.out.getByteCount () + " bytes";
}, "~S");
});
