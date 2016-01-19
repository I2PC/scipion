// Scipion JS app: this must be populated by the templates: always!!
var scipion = {
    'subdomain': ''
}

function setSubdomainURL(subdomain){

    // Everycode that uses getSubdomainURL() expect NOT to have a backslash..so remove it.
    // Get the last character
    var last = subdomain.slice(-1);

    if (last == "/") {
        subdomain = subdomain.slice(0, -1);
    }

    scipion.subdomain = subdomain;

}

function getSubDomainURL(){return scipion.subdomain}

function getPort(){
	if (document.location.port != "") {
		return ":" + document.location.port;
	} else {
		return "";
	}
}
function getBaseURL() {

	return "http://" + document.domain + getPort() + getSubDomainURL();
}

function getAbsoluteURL(url){
	var domain = getSubDomainURL()
	if (url.indexOf(domain + "/") != 0) {
		var url = domain + url
    }
	return url
}

function goWithSubDomainURL(url){window.location.href=getAbsoluteURL(url)}
function getFormUrl(){return $("input#formUrl").val()}
function setFormURL(url){
	var formURL = getFormUrl();
	var newURL = url.split("/form/")
	return "/"+formURL+"/"+newURL[1]
}