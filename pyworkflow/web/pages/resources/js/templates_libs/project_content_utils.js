/**
 * Methods used in the project_content template. 
 * Toolbar functions + Manage tabs
 * 
 * launchToolbar(projName, id, elm);
 * fillTabsSummary(id);
 * fillUL(list, ulId, icon);
 * launchViewer(id);
 * updateTabs(projName, id, elm);
 * switchGraph();
 * deleteProtocolForm(projName, protocolId);
 * deleteProtocol(elm);
 *  
 * 
 **/

/*
 * Toolbar used in the project content template for list view
 */
function launchToolbarList(projName, id, elm) {
	var row = $("div#toolbar");
	updateRow(id, elm, row);
	updateButtons(projName, id, elm);
	row.show(); // Show toolbar
}

/*
 * Toolbar used in the project content template for list view
 */
function launchToolbarTree(projName, id, elm) {
	var row = $("div#toolbar");
	updateTree(id,elm);
	updateButtons(projName, id, elm);
	row.show(); // Show toolbar
}

function checkRunStatus(projName, id) {
	$.ajax({
		type : "GET",
		url : '/protocol_status/?protocolId=' + id,
		dataType:"text",
		success : function(status) {
			if(status=="running"){
				// Action Stop Button
				$("span#analyzeTool").hide();
				$("span#stopTool").show();
				$("a#stopTool").attr('href',
				'javascript:stopProtocolForm("' + projName + '","' + id + '")');
			}
			else{
				// Action Analyze Result Button
				$("span#stopTool").hide();
				$("span#analyzeTool").show();
				$("a#analyzeTool").attr('href', 'javascript:launchViewer("'+id +'")');
			}
		}
	});
}

/*
 * Fill the content of the summary tab
 */
function fillTabsSummary(id) {
	$.ajax({
		type : "GET",
		url : '/protocol_io/?protocolId=' + id,
		dataType : "json",
		success : function(json) {
			fillUL(json.inputs, "protocol_input", "db_input.gif");
			fillUL(json.outputs, "protocol_output", "db_output.gif");
		}
	});

	$.ajax({
		type : "GET",
		url : '/protocol_summary/?protocolId='+ id,
		dataType : "json",
		success : function(json) {
			$("#tab-summary").empty();
			for ( var i = 0; i < json.length; i++) {
				$("#tab-summary").append('<p>' + json[i] + '</p>');
			}
		}
	});
}

/*
 * Fill an UL element with items from a list items should contains id and name
 * properties
 */
function fillUL(list, ulId, icon) {
	ul = $("#" + ulId);
	ul.empty();
	for ( var i = 0; i < list.length; i++) {
		ul.append('<li><a href="/visualize_object/?objectId=' + list[i].id
				+ '"target="_blank"><img src="../../../../resources/' + icon + '" /> '
				+ list[i].name + '</a></li>');
	}
}

/*
 * Launch the viewers to analyze the results of the protocol run
 */
function launchViewer(id){
	/* Execute the viewer */
	$.ajax({
		type : "GET",
		url : "/viewer/?protocolId=" + id ,
		dataType : "json",
		success : function(json) {
			$.each(json, function(key, value) {
				if(value.substring(0,6)=="<html>"){
					//Launch a pop up with the html
					customPopupHTML(value,800,600);
				}else{
					//Launch a pop up with the url
					customPopup(value,1024,768);
				}
			});
		}
	});	
}

function updateButtons(projName, id, elm){
	// Action Edit Button
	$("a#editTool").attr('href',
	'javascript:popup("/form/?=&protocolId=' + id + '")');
	
	// Action Copy Button
	$("a#copyTool").attr('href',
	'javascript:popup("/form/?&protocolId=' + id + '&action=copy' + '")');

	// Action Delete Button
	$("a#deleteTool").attr('href',
			'javascript:deleteProtocolForm("' + projName + '","' + id + '")');

	// Action Browse Button
	var aux = "javascript:alert('Not implemented yet')";
	$("a#browseTool").attr('href', aux);
	
	checkRunStatus(projName, id);
	fillTabsSummary(id);
}


function updateTree(id, elm){
	var oldSelect = $("div#graphActiv").attr("data-option");
	var selected = "graph_" + id;

	if (oldSelect != selected) {
		if (oldSelect != "") {
			var aux = "div#" + oldSelect + ".window";
			$(aux).css("border", "");
		}
		$("div#graphActiv").attr("data-option", selected);
		elm.css("border", "2.5px solid Firebrick");
	}
}
	
function updateRow(id, elm, row){	
	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', 'runtr');
	}
	row.attr('value', id);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');
}


function switchGraph() {
	var status = $("div#graphActiv").attr("data-mode");
	// Graph will be painted once
	if ($("div#graphActiv").attr("data-time") == 'first') {
		if (status == 'inactive') {
			// Graph ON
			$("div#graphActiv").attr("data-mode", "active");
			$("div#graphActiv").attr("style", "");
			$("div#treeTool").hide();
			// Table OFF
			$("div#runTable").attr("data-mode", "inactive");
			$("div#runTable").attr("style", "display:none;");
			$("div#listTool").show();
		} else if (status == 'active') {
			// Table ON	
			$("div#runTable").attr("data-mode", "active");
			$("div#runTable").attr("style", "");
			$("div#listTool").hide();
			// Graph OFF
			$("div#graphActiv").attr("data-mode", "inactive");
			$("div#graphActiv").attr("style", "display:none;");
			$("div#treeTool").show();
		}
		callPaintGraph();
		$("div#graphActiv").attr("data-time", "not");
		
	} else {
		if (status == 'inactive') {
			// Graph ON
			$("div#graphActiv").attr("data-mode", "active");
			$("div#graphActiv").attr("style", "");
			$("div#treeTool").hide();
			// Table OFF
			$("div#runTable").attr("data-mode", "inactive");
			$("div#runTable").attr("style", "display:none;");
			$("div#listTool").show();
			
			// getElement in table
			var s = $("tr.selected").attr("id");
			s = "graph_" + s;

			if (s != "") {
				var nodeClear = $("div#graphActiv").attr("data-option");
				if (nodeClear != "") {
					if (nodeClear != s) {
						// Clear the node selected
						var elmClear = $("div#" + nodeClear + ".window");
						elmClear.css("border", "");

						// setElement in graph
						$("div#graphActiv").attr("data-option", s);

						// Highlight the node
						var elm = $("div#" + s + ".window");
						elm.css("border", "2.5px solid Firebrick");

					}
				}
			}

		} else if (status == 'active') {
			// Table ON
			$("div#runTable").attr("data-mode", "active");
			$("div#runTable").attr("style", "");
			$("div#listTool").hide();
			// Graph OFF
			$("div#graphActiv").attr("data-mode", "inactive");
			$("div#graphActiv").attr("style", "display:none;");
			$("div#treeTool").show();
			
			// getElement in graph
			var s = $("div#graphActiv").attr("data-option");
			var s = s.replace("graph_", "");

			if (s != "") {
				var rowClear = $("tr.selected").attr("id");
				if (rowClear != "") {
					if (rowClear != s) {
						// Clear the row selected
						var elmClear = $("tr.selected");
						elmClear.attr("style", "background-color: #fafafa;");
						elmClear.attr("class", "runtr");

						// setElement in table
						var elm = $("tr#" + s + ".runtr");
						var projName = $("div#graphActiv").attr("data-project");
						// elm.attr("style", "background-color:
						// LightSteelBlue;");
						// elm.attr("class","selected");
						launchToolbarList(projName, s, elm);
					}
				}
			}
		}
	}
}

/*
 * Dialog form to verify the right option to delete
 */
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

/*
 * Method to execute a delete by a protocol
 */
function deleteProtocol(elm) {
	var value = elm.attr('value').split("-");
	var projName = value[0];
	var protId = value[1];
	$.ajax({
		type : "GET",
		url : "/delete_protocol/?protocolId=" + protId
	});
}

/*
 * Dialog form to verify the right option to stop a protocol
 */
function stopProtocolForm(projName, protocolId) {

	var msg = "<table><tr><td><img src='/resources/warning.gif' width='45' height='45' />"
			+ "</td><td class='content' value='"
			+ projName
			+ "-"
			+ protocolId
			+ "'>This <strong>protocol run</strong>"
			+ " will be <strong>STOPPED</strong>. Do you really want to continue?</td></tr></table>";

	new Messi(msg, {
		title : 'Confirm STOP',
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'btn-select',
			btnFunc : 'stopProtocol'
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

/*
 * Method to stop the run for a protocol
 */
function stopProtocol(elm) {
	var value = elm.attr('value').split("-");
	var projName = value[0];
	var protId = value[1];
	$.ajax({
		type : "GET",
		url : "/stop_protocol/?protocolId=" + protId
	});
}
