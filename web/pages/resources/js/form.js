function evalElements() {
	var elements = document.getElementsByTagName("tr");
	var size = elements.length;
	
	alert(size);
	

	for ( var x = 0; x < size; x++) {
		
		alert(elements.item(x).attr('name'));
//		 alert(elements.item(x).textContent);
//		 evalDependencies(elements.item(x));
	}
}

function evalDependencies(aux) {
	alert(aux);
	var name = aux.attr('name');
	var row = jQuery("tr#" + name);
	var dependencies = row.attr('data-depen');
	if (dependencies.length > 0) {
		var arrayDepends = dependencies.split(",");
		for ( var cont = 0; cont < arrayDepends.length; cont++) {
			var res = evalCondition(arrayDepends[cont]);
			if (res == false) {
				jQuery("tr#" + arrayDepends[cont]).hide();
			} else if (res == true) {
				jQuery("tr#" + arrayDepends[cont]).show();
			}
		}
	}
}

function getValueByName(itemName) {
	var row = jQuery("tr#" + itemName);
	/*
	 * Get the type of the input with the next cases: BooleanParam, StringParam,
	 * EnumParam (type select or inputs)
	 */
	var type = row.attr('data-type');
	var value = null;

	if (type == 'BooleanParam') {
		value = jQuery("input#" + itemName + "_yes").attr('checked')
	} else if (type == 'StringParam') {
		value = jQuery("input#" + itemName + "_input").attr('value');
	} else if (type == 'FloatParam') {
		value = jQuery("input#" + itemName + "_input").attr('value');
	} else if (type == 'EnumParam') {
		var enumType = row.attr('data-enum');
		if (enumType == '1') {
			var elm = jQuery("select#" + itemName + "_select");
			value = elm.attr('value');
		} else if (enumType == '0') {
			var cont = 0;
			var enc = 0;
			while (!enc) {
				var opt = jQuery("input#" + itemName + "_" + cont);
				if (opt.attr('checked') == true) {
					value = cont;
					value = opt.attr('value');
					enc = 1;
				}
				cont++;
			}
		}
	}
	// alert("getValue: " + itemName + "=" + value);
	return value;
}

function evalCondition(itemName) {
	var row = jQuery("tr#" + itemName);
	var type = row.attr('data-type');
	var cond = row.attr('data-cond');
	var params = row.attr('data-params');

	var arrayParams = params.split(",");

	// Get value of the element with name=itenName
	var param = null;
	var value = null;
	var cond_eval = cond;

	for ( var cont = 0; cont < arrayParams.length; cont++) {
		param = arrayParams[cont];
		value = getValueByName(param);
		cond_eval = cond_eval.replace(param, value);
	}

	// alert("condition: " + cond + " eval: " + cond_eval);
	return eval(cond_eval);
}

function help(title, msg) {
	new Messi(msg, {
		title : 'Help' + ' ' + title,
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Close',
			val : 'X',
			btnClass : 'btn-close'
		} ]
	});
}

function browse(title, msg) {
	new Messi('', {
		title : msg,
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Choose',
			val : 'Y',
			btnClass : 'btn-choose'
		}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'btn-cancel'
		} ],
		callback : function(val) {
			if (val == 'Y') {
				alert('You are selected one');
			}
		}
	});
}

function filemanager(elm, type, name) {
	alert(elm + " " + type + " " + name);
	// document.getElementById(elm).innerText = document
	// .getElementById('openssme').value
}
