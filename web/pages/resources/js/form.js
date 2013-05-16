function evalCondition(aux) {
	var name = aux.attr('name')
	var value = aux.attr('value')
	var affected = aux.attr('data-ref')
	var array = affected.split(",");

	if (value == 'no') {
		for (cont = 0; cont < array.length; cont++) {
			jQuery("tr#" + array[cont]).hide();
		}
	} else if (value == 'yes') {
		for (cont = 0; cont < array.length; cont++) {
			jQuery("tr#" + array[cont]).show();
		}
	}
	// alert(b.getAttribute("data-ref"));
}

// function evalCondition(aux1, aux2) {
// var name = aux.attr('name')
// var value = aux.attr('value')
// var affected = aux.attr('data-ref')
//
// if (value == 'no') {
// jQuery("tr#" + affected).hide();
// } else if (value == 'yes') {
// jQuery("tr#" + affected).show();
// }
// // alert(b.getAttribute("data-ref"));
// }
//
// function evalExpLevel(aux1, aux2) {
// var name = aux.attr('name')
// var value = aux.attr('value')
// var affected = aux.attr('data-ref')
//
// if (value == 'no') {
// jQuery("tr#" + affected).hide();
// } else if (value == 'yes') {
// jQuery("tr#" + affected).show();
// }
// // alert(b.getAttribute("data-ref"));
// }
