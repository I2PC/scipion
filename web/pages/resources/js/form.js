function evalCondition(aux) {
	var name = aux.attr('name')
	var value = aux.attr('value')
	var show = aux.attr('data-show')
	var hide = aux.attr('data-hide')

	if (show != null) {
		var array = show.split(",");

		for (cont = 0; cont < array.length; cont++) {
			jQuery("tr#" + array[cont]).show();
		}
	}

	if (hide != null) {
		var array = hide.split(",");

		for (cont = 0; cont < array.length; cont++) {
			jQuery("tr#" + array[cont]).hide();
		}
	}
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
