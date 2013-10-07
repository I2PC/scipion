/**
 * Methods to manage the protocol form
 * 
 * showErrorValidation(json);
 * evalElements();
 * onChangeParam(value, paramId);
 * onChangeEnumParamCombo(elemId, paramId);
 * onChangeEnumParamList(index, paramId);
 * setParamValue(paramId, value);
 * evalDependencies(row, newLevel);
 * evalCondition(row);
 * help(title, msg);
 * browseObjects(node, projName, objClass);
 * getListFormatted(node, list, id);
 * getTableFormatted(node, list, id);
 * selectDialog(objClass, msg, funcName);
 * processSelectionList(elm);
 * processSelectionTable(elm);
 * selTableMessi(elm);
 * 
 **/

$(document).ready(function() {
	/*	
		 * Method to execute a protocol.
		 * Overray the post simple method in the html. 
	*/
	$("#protocolForm").submit(function() {
		var mode = $("#protocolForm").attr('data-mode');

		if (mode == 'execute') {
			/* Execute the protocol */
			var action = $("#protocolForm").attr("action");
			
			var msg = messiInfo("The protocol was launched successfuly");

			$.post(action, $("#protocolForm").serialize(), function(json) {
				if (json.errors.length > 0) {
					// Show errors in the validation
					showErrorValidation(json.errors);
				} else {
					// No errors in the validation
					new Messi(msg, {
						title : 'Success',
						modal : true,
						buttons : [ {
							id : 0,
							label : 'Ok',
							val : 'Y',
							btnClass : 'btn-select'
						} ],
						callback : function(val) {
							if (val == 'Y') {
								window.opener.location.reload(true);
								window.close();
							}
						}
					});
				}
			}, "json");

		} else if (mode == 'save') {
			/* Save the protocol */
			var action = "/save_protocol/";
			
			var msg = messiInfo("The protocol was saved successfuly");

			$.post(action, $("#protocolForm").serialize(), function() {
				new Messi(msg, {
					title : 'Success',
					modal : true,
					buttons : [ {
						id : 0,
						label : 'Ok',
						val : 'Y',
						btnClass : 'btn-select'
					} ],
					callback : function(val) {
						if (val == 'Y') {
							window.opener.location.reload(true);
							// window.close();
						}
					}
				});
			});
		} else if (mode == 'wiz') {
			
			new Messi("<img src='/resources/tools_wizard.png'/>  Loading Wizard...",{
				modal : true
				});
			
			/* Execute the wizard */
			var action = "/wizard/";
			var type_wiz = $("#wizName").attr("value");
			
			$.post(action, $("#protocolForm").serialize(), function(html) {
				
				$('.messi').remove();
				$('.messi-modal').remove();
				
				if(html=="errorInput"){
					var msg = messiError("Input was not selected, please choose one.");
					launchMessiSimple("Error",msg);
				} else if (html=="errorEmpty"){
					var msg = messiError("Input particles selected are None");
					launchMessiSimple("Error",msg);
				} else if (html=="errorIterate"){
					var msg = messiError("Error iterating over the set of particles");
					launchMessiSimple("Error",msg);
				} else if(type_wiz=='wiz_particle_mask' || type_wiz=='wiz_volume_mask'){
					customPopupHTML(html,520,460);
				} else if(type_wiz=='wiz_volume_mask_radii'){
					customPopupHTML(html,530,520);
				} else{
					customPopupHTML(html,775,470);
				}
			});
		}
		// Important. Stop the normal POST
		return false;
	});
});

function showErrorValidation(json) {
	var msg = JSON.stringify(json);
	msg = msg.replace("<", "");
	msg = msg.replace(">", "");
	msg = msg.replace("[", "");
	msg = msg.replace("]", "");

	var msg = "<table><tr><td><img src='/resources/error.gif' width='45' height='45' />"
			+ "</td><td class='content'>" + msg + "</td></tr></table>";

	new Messi(msg, {
		title : 'Errors found',
		modal : true,
		buttons : [ {
			id : 0,
			label : 'close',
			val : 'Y',
			btnClass : 'btn-close'
		} ]
	});
}

function evalElements() {
	$("tr").each(function(index) {
//		
		var value = jQuery(this).val();
		if(value.length == 0){
			var value = jQuery(this).attr('value');
		}
		var type = jQuery(this).attr('data-type');
		var param = jQuery(this).attr('id');

//		 alert(value +" - "+type+" - "+param);

//		if (type == "BooleanParam" || type == "FloatParam" || type == "IntParam") {
//			onChangeBooleanParam(value, param);
//		} else 
		if (type == "EnumParam") {
			var typeEnum = jQuery(this).attr('data-enum');
			if (typeEnum == '0') {TATOO
				onChangeEnumParamList(value, param);
			} else if (typeEnum == '1') {
				onChangeEnumParamCombo(param + "_select", param);
			}
		} else {
//			alert("before:"+value);
			onChangeParam(value, param);
		}
	});
}

/* Differents functions depends on the input type */
function onChangeParam(value, paramId) {
//	alert(paramId + "-"+value);
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
//	alert("before:"+row.val());
	row.val(value);	
//	alert("after:"+row.val());
	
	var newLevel = $("select[name=expertLevel]").val();
	evalDependencies(row, newLevel);

	var params = row.attr('data-params');

	if (params != undefined && params.length <= 0) {
					
		var expLevel = row.attr('data-expert');
	
		if (expLevel > newLevel) {
			row.hide();
		} else {
			row.show();			
		}
	}
}

function evalDependencies(row, newLevel) {
//	var newLevel = $("select[name=expLevel]").val();

	var dependencies = row.attr('data-depen');
	if (dependencies!= undefined && dependencies.length > 0) {
		var arrayDepends = dependencies.split(",");
		for ( var cont = 0; cont < arrayDepends.length; cont++) {

			var row2 = jQuery("tr#" + arrayDepends[cont]);
			var res = evalCondition(row2);
			var expLevel = row2.attr('data-expert');

			if (res == false || expLevel > newLevel) {
				row2.hide();
			} else if (res == true) {
				row2.show();
				evalDependencies(row2, newLevel);
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
		if (!value){
			value = jQuery("tr#" + param).attr("value");
			if (!value){
				value="''";
			}
		}
		params += "param: " + param + " value: " + value + "\n";
		cond_eval = cond_eval.replace(param, value);
	}
	// if (row.attr("name")=="comment") {
	// alert("condition: " + cond + " \nparams:\n" + params + "\n eval: " +
	// cond_eval);
	// }
	
	cond_eval = normalizeConditions(cond_eval)

	var foundAnd = cond.indexOf("'0'") != -1;
	if (foundAnd)
		alert("condition: " + cond + " \neval: " + cond_eval+ " \nparams:"+ params );

	return eval(cond_eval);
}

function normalizeConditions(cond){
	cond = cond.replace("not","!");
	cond = cond.replace("and","&&");
	cond = cond.replace("or","||");
//	cond = cond.replace("'0'","false");
//	cond = cond.replace("'1'","true");
	return cond;
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

function showComment() {
	var msg = $("input#comment").attr("value");
	
	if(msg == ""){
		msg = "Describe your run here...";
	}
	
	var msg ="<textarea id='description'>"+ msg +"</textarea>";
	
	new Messi(msg, {
		title : 'Comment',
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Select',
			val : 'Y',
			btnClass : 'btn-select',
			btnFunc : 'putComment'
			}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'btn-cancel'
		}]
	});
}

function putComment(){
	$("input#comment").attr("value",$("textarea#description").val());
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
		title : 'Select ' + objClass,
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

/*
 * Used to choose a element in the protocol form
 */
function selTableMessi(elm) {

	var row = $("table.content");
	var id = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("td#" + row.attr('value'));
		rowOld.attr('style', '');
	}
	row.attr('value', id);
	elm.attr('style', 'background-color: LightSteelBlue;');
}

