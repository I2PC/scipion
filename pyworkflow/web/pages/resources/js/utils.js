function popup(URL) {
	var popup_width = 470
	var popup_height = 615
	day = new Date();
	id = day.getTime();
	eval("page"
			+ id
			+ " = window.open(URL, '"
			+ id
			+ "', 'toolbar=0,scrollbars=1,location=0,statusbar=0,menubar=0,resizable=0,width='+popup_width+',height='+popup_height+'');");
}

function customPopup(URL, widthValue, heightValue) {
	day = new Date();
	id = day.getTime();
	eval("page"
			+ id
			+ " = window.open(URL, '"
			+ id
			+ "', 'toolbar=0,scrollbars=1,location=0,statusbar=0,menubar=0,resizable=0,width='+widthValue+',height='+heightValue+'');");
}

function closePopup() {
//	opener.location.reload(true);
//	self.close();
	window.opener.location.reload(true);
	window.close();
}

/*
 * Toolbar used in the project content template
 */
function launchToolbar(projName, id, elm) {
	var row = $("div#toolbar");

	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', '');
	}
	row.attr('value', id);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');

	// Action Edit Button
	$("a#editTool").attr(
			'href',
			'javascript:popup("/form/?projectName=' + projName + '&protocolId='
					+ id + '")');
	// Action Copy Button
	$("a#copyTool").attr(
			'href',
			'javascript:popup("/form/?projectName=' + projName + '&protocolId='
					+ id + '&action=copy' + '")');
	// Action Delete Button
	$("a#deleteTool").attr('href',
			'javascript:deleteProtocolForm("' + projName + '","' + id + '")');
	// Action Browse Button
	// $("a#browseTool").attr(
	// 'href',
	// 'javascript:popup("/form/?projectName=' + projName + '&protocolId='
	// + id + '")');
	//alert('launchToolbar')
	
	row.show(); // Show toolbar
	
	$.ajax({
		type : "GET",
		url : '/protocol_io/?projectName=' + projName + '&protocolId='
					+ id,
		dataType : "json",
		success : function(json) {
			fillUL(json.inputs, "protocol_input", "db_input.gif", projName)
			fillUL(json.outputs, "protocol_output", "db_output.gif", projName)
			
//			ul_output = $("#protocol_output")
//			ul_output.empty()
//			for (var i = 0; i < json.outputs.length; i++) {
//				ul_output.append(
//				'<li><a href=""><img src="../../../../resources/db_output.gif" /> ' + json.outputs[i].name
//				+ '</a></li>');				
//			}
			//'<li><a href="/user/messages"><span class="tab">Message Center</span></a></li>');
//				for ( var x = 0; x < list.length; x++) {
//					res += "<input type='radio' id ='" + id + x + "' name='" + id
//							+ "'  value='" + list[x] + "' />" + list[x] + "<br />";
//				}
		}
	});
}

/*
 * Fill an UL element with items from a list
 * items should contains id and name properties
 */
function fillUL(list, ulId, icon, projName) {
	ul = $("#" + ulId)
	ul.empty()
	for (var i = 0; i < list.length; i++) {
		ul.append(
				'<li><a href="/visualize_object/?projectName=' + projName + 
				'&objectId=' + list[i].id +
				'"><img src="../../../../resources/' + icon + '" /> ' + list[i].name
				+ '</a></li>');				
	}		
}

/*
 * Toolbar used in the view host template
 */
function launchHostsToolbar(projName, hostId, elm) {
	var row = $("div#toolbarHost");

	if (row.attr('value') != undefined && row.attr('value') != hostId) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', '');
	}
	row.attr('value', hostId);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');

	// Action Edit Button
	$("a#editTool").attr(
			'href',
			'javascript:editHost()');
	// Action Copy Button
	$("a#newTool").attr(
			'href',
			'javascript:newHost()');
	// Action Delete Button
	$("a#deleteTool").attr('href',
			'javascript:deleteHost()');
	// Action Browse Button
	// $("a#browseTool").attr(
	// 'href',
	// 'javascript:popup("/form/?projectName=' + projName + '&protocolId='
	// + id + '")');

	row.show(); // Show toolbar
}

function deleteProtocolForm(projName, protocolId) {

	var msg = "<table><tr><td><img src='/resources/warning.gif' width='45' height='45' />"
			+ "</td><td class='content' value='"
			+ projName
			+ "-"
			+ protocolId
			+ "'><strong>ALL DATA</strong> related to this <strong>protocol run</strong>"
			+ " will be <strong>DELETED</strong>. Do you really want to continue?</td></tr></table>";

	new Messi(msg, {
		title : 'Confirm DELETE',
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'btn-select',
			btnFunc : 'deleteProtocol'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'btn-cancel'
		} ],
		callback : function(val) {
			if (val == 'Y') {
				window.location.href = "/project_content/?projectName="
						+ projName;
			}
		}
	});
}

function deleteProtocol(elm) {
	var value = elm.attr('value').split("-");
	var projName = value[0];
	var protId = value[1];
	$.ajax({
		type : "GET",
		url : "/delete_protocol/?projectName=" + projName + "&protocolId="
				+ protId
	});
}

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
