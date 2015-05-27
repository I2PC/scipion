function getSubDomainURL(){return ""}
function setSubDomainURL(url){
	var domain = getSubDomainURL() + "/"
	if (url.indexOf(domain) != 0) {
		var url = domain + url
    }
	return url
}
function goWithSubDomainURL(url){window.location.href=setSubDomainURL(url)}
function getFormUrl(){return $("input#formUrl").val()}
function setFormURL(url){
	var formURL = getFormUrl();
	var newURL = url.split("/form/")
	return "/"+formURL+"/"+newURL[1]
}