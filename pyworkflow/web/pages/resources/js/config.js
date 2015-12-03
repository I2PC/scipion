function getSubDomainURL(){return ""}

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