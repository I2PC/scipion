function evalElements() {
	alert("hola");
	$("tr").each(function(index) {
		var value = jQuery(this).attr('data-value');
		var type = jQuery(this).attr('data-type');
		var param = jQuery(this).attr('id');

		if (type == "BooleanParam") {
			onChangeBooleanParam(value, param);
		} else if (type == "EnumParam") {
			var typeEnum = jQuery(this).attr('data-enum');
			if (typeEnum == '0') {
				onChangeEnumParamList(value, param);
			} else if (typeEnum == '1') {
				onChangeEnumParamCombo(param + "_select", param);
			}
		}
	});
}

/* Differents functions depends on the input type */
function onChangeBooleanParam(value, paramId) {
	setParamValue(paramId, value);
}

function onChangeEnumParamCombo(elemId, paramId) {
	var elem = document.getElementById(elemId);
	setParamValue(paramId, elem.selectedIndex);
}

function onChangeEnumParamList(index, paramId) {
	setParamValue(paramId, index);
}

// Put the new value in an attribute of the parent node
function setParamValue(paramId, value) {
	var row = jQuery("tr#" + paramId);
	row.val(value);
	evalDependencies(row);
}

function evalDependencies(row) {
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

function evalCondition(itemName) {
	var row = jQuery("tr#" + itemName);
	var cond = row.attr('data-cond');
	var params = row.attr('data-params');

	var arrayParams = params.split(",");

	// Get value of the element with name=itenName
	var param = null;
	var value = null;
	var cond_eval = cond;

	for ( var cont = 0; cont < arrayParams.length; cont++) {
		param = arrayParams[cont];
		// value = getValueByName(param);
		value = jQuery("tr#" + param).val();
		cond_eval = cond_eval.replace(param, value);
	}
	// alert("condition: " + cond + " eval: " + cond_eval);
	if (cond_eval == "True") {
		return true;
	} else if (cond_eval == "False") {
		return false;
	} else {
		return eval(cond_eval);
	}
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

/*
 * Browse object in the database. Params: objClass: the class to get instances
 * from (also subclasses)
 */
function browseObjects(objClass) {
	$.ajax({
		type : "GET",
		url : "/browse_objects/?objClass=" + objClass,
		dataType : "json",
		success : function(json) {
			// specifying a dataType of json makes jQuery pre-eval the response
			// for us
			// alert(json.objects);

			//Generate the list formatted
			var array = json.objects;
			var res = getTableFormatted(array);

			selectObjects('Select ' + objClass, res);

		}
	});
}

function getTableFormatted(list){
	var res = "<div style='overflow:auto'>
	<table class='browse' id='browse' cellspacing='0'>";
	for(var x=0;x<list.length;x++){
		res = res + "<tr id='browse'>
		<td id='browse'><a id='browse' href=''>"+ list[x] + 
		"</a></td></tr>";
	}
	res = res + "</table></div>";
	return res;
}

function selectObjects(title, msg) {
	new Messi(msg, {
		title : title,
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Select',
			val : 'Y',
			btnClass : 'btn-select'
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
