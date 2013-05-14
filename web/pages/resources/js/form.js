
function evalCondition(aux) {
	var name = aux.attr('name')
	var value = aux.attr('value')
	var affected = aux.attr('data-ref')

//	var b = document.getElementById(affected);

	if (value = 'no') {
		// b.setAttribute("style", "display:none;");
		// b.style.visibility = "hidden";
		jQuery("div#" + affected).hide();
	} else if (value = 'yes') {
		// b.setAttribute("style", "display:block;")
		// b.style.visibility = "visible";
		jQuery("div#" + affected).show();
	}

	// 
	// alert(b.getAttribute("data-ref"));
	alert(name + " " + value);
}