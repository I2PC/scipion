function evalElements() {
	$("tr").each(function(index) {
		var value = jQuery(this).attr('value');
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
	var newLevel = $("select[name=expLevel]").val();

	var dependencies = row.attr('data-depen');
	if (dependencies.length > 0) {
		var arrayDepends = dependencies.split(",");
		for ( var cont = 0; cont < arrayDepends.length; cont++) {

			var row = jQuery("tr#" + arrayDepends[cont]);
			var res = evalCondition(row);
			var expLevel = row.attr('data-expert');

			if (res == false || expLevel > newLevel) {
				row.hide();
			} else if (res == true) {
				row.show();
			}
		}
	}
}

function evalCondition(row) {

	var cond = row.attr('data-cond');
	var params = row.attr('data-params');

	var arrayParams = params.split(",");

	// Get value of the element with name=itenName
	var param = null;
	var value = null;
	var cond_eval = cond;
	// var params = '';

	for ( var cont = 0; cont < arrayParams.length; cont++) {
		param = arrayParams[cont];
		// value = getValueByName(param);
		value = jQuery("tr#" + param).val();
		// params += "param: " + param + " value: " + value + "\n";
		cond_eval = cond_eval.replace(param, value);
	}
	// if (row.attr("name")=="comment") {
	// alert("condition: " + cond + " \nparams:\n" + params + "\n eval: " +
	// cond_eval);
	// }

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

/*
 * Browse object in the database. Params: objClass: the class to get instances
 * from (also subclasses)
 */
function browseObjects(node, projName, objClass) {
	$.ajax({
		type : "GET",
		url : "/browse_objects/?projectName=" + projName + "&objClass="
				+ objClass,
		dataType : "json",
		success : function(json) {
			// specifying a dataType of json makes jQuery pre-eval the response
			// for us
			// var res = getListFormatted(node, json.objects, objClass);
			var res = getTableFormatted(node, json.objects, objClass);

			// selectDialog(objClass, res, "processSelectionList");
			selectDialog(objClass, res, "processSelectionTable");

		}
	});
}

function getListFormatted(node, list, id) {
	var res = "<div class='content' style='overflow:auto' data-node='" + node
			+ "'>";
	for ( var x = 0; x < list.length; x++) {
		res += "<input type='radio' id ='" + id + x + "' name='" + id
				+ "'  value='" + list[x] + "' />" + list[x] + "<br />";
	}
	res = res + "</div>";
	return res;
}

function getTableFormatted(node, list, id) {

	var res = "<table class='content' style='overflow:auto' data-node='" + node
			+ "'>";
	for ( var x = 0; x < list.length; x++) {
		res += "<tr><td id='" + id + x + "' name='" + id + "' value='"
				+ list[x] + "' onclick=javascript:selTableMessi($(this)); >"
				+ list[x] + "</td></tr>";
	}
	res = res + "</table>";
	return res;
}

function selectDialog(objClass, msg, funcName) {
	new Messi(msg, {
		title : 'Select' + objClass,
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Select',
			val : 'Y',
			btnClass : 'btn-select',
			btnFunc : funcName
		}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'btn-cancel'
		} ]
	// callback : function(val) {
	// if (val == 'Y') {
	// alert();
	// }
	// }
	});
}

function processSelectionList(elm) {
	elm.children('input').each(function() {
		if (jQuery(this).attr('checked')) {
			elm.val(jQuery(this).attr('value'));
		}
	});
	jQuery('input#' + elm.attr('data-node') + '_input').val(elm.attr('value'));
}

function processSelectionTable(elm) {
	var selected = elm.attr('value');
	var value = $("td#" + selected).attr('value');

	jQuery('input#' + elm.attr('data-node') + '_input').val(value);
}
