
function evalCondition(aux) {
	var name = aux.attr('name')
	var value = aux.attr('value')
	var affected = aux.attr('data-ref')

	var b = document.getElementById(affected);

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
	// alert(name + " " + value);

}

function invisibleRapida() {
	var a = document.getElementById("busquedaRapida");
	a.style.visibility = "hidden";
}

function invisibleAvanzada() {
	var b = document.getElementById("busquedaAvanzada");
	b.setAttribute("style", "display:none;");
	b.style.visibility = "hidden";
}

function visibleRapida() {
	var a = document.getElementById("busquedaRapida");
	a.style.visibility = "visible";
}

function visibleAvanzada() {
	var b = document.getElementById("busquedaAvanzada");
	b.setAttribute("style", "display:block;")
	b.style.visibility = "visible";
}

function cambioRapAv() {
	invisibleRapida();
	visibleAvanzada();
}

function cambioAvRap() {
	invisibleAvanzada();
	visibleRapida();
}

function visibleRapidaHTML() {
	jQuery('div.form_avanzada').hide();
	jQuery('div.form_rapida').show();

}

function visibleAvanzadaHTML() {
	jQuery('div.form_rapida').hide();
	jQuery('div.form_avanzada').show();
}

function cambioAvRapHTML() {
	visibleRapidaHTML();
}

function cambioRapAvHTML() {
	visibleAvanzadaHTML();
}

function invisibleResultado() {
	jQuery('div.resultado').hide();
}
