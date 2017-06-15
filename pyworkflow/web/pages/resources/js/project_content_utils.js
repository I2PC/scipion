/*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'scipion@cnb.csic.es'
 *
 ******************************************************************************/
/******************************************************************************
 * DESCRIPTION:
 * 
 * Methods used in the project_content template. 
 * Toolbar functions + Manage tabs
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * 
 * function preloadToolbar(list)
 * 	->	Function to preload the nodes/rows selected.
 * 
 * function launchToolbarList(id, elm)
 * 	->	Toolbar used in the project content template for the run list view.
 * 
 * function launchToolbarTree(id, elm, multi)
 * 	->	Toolbar used in the project content template for graph view.
 * 
 * function launchToolbarProject(id, elm, type)
 * 	->	Method used to update an element in the template depending on 
 * 		view mode (list or graph).
 * 
 * function refreshToolbarSingleMark(id)
 * 	->	Update the buttons functionalities into the toolbar
 * 		for single marks.
 * 
 * function refreshToolbarMultipleMark(id, elm)
 * 	->	Update the buttons functionalities into the toolbar
 * 		for multiple marks.
 * 
 * function enableMultipleMarkGraph(elm)
 * 	->	Used to highlight a node in the protocol graph.
 * 		This method is used to multiples marked nodes.
 * 
 * function disableMultipleMarkGraph(elm)
 * 	->	Used to remove the highlight applied to a node in the protocol graph.
 * 		This method is used to multiples marked nodes.
 * 
 * function getElmMarkedList()
 * 	->	Return a list of elements marked in the view mode currently.
 * 
 * function getElmMarkedGraph()
 * 	->	Return a list of elements marked in the graph view mode.
 * 
 * function transposeElmMarked(status)
 * 	->	Switch the elements marked between graph and list view.
 * 
 * function refreshSelectedRuns(list_marked)
 * 	->	Refresh the elements contained into the list.
 * 
 * --- GRAPH METHODS ---
 * function markSingleNodeGraph(elm)
 * ->	Method to mark a node selected in the graph.
 * 
 * dismarkSingleNodeGraph(elm)
 * 	->	Method to dismark a node selected in the graph.
 * 
 * function enableMultipleMarkGraph(elm)
 * 	->	Method to active the multiple selection.
 * 
 * function disableMultipleMarkGraph(id)
 * 	->	Method to disactive the multiple selection.
 * 
 * --- LIST METHODS --- 
 * function markSingleNodeList(elm)
 * 	->	Method to mark a row selected in the list.
 * 
 * function dismarkSingleNodeList(elm)
 * 	->	Method to dismark a row selected in the list.
 * 
 * function enableMultipleMarkList(elm)
 * 	->	Method to activate the multiple selection.
 * 
 * function disableMultipleMarkList(id)
 * 	->	Method to disactivate the multiple selection.
 * --------------------
 * 
 * function checkRunStatus(id)
 * 	->	Function to check a protocol run, depend on the status two button will be
 * 		switching in the toolbar (Stop / Analyze Results).
 * 
 * function updateTabs(id)
 * 	->	Fill the content of the tabs for a protocol run selected 
 * 		(Data / Summary / Methods / Status)
 * 
 * function showLog(log_type)
 * 	->	This function is used to show or hide differents logs about a protocol selected.
 * 
 *  * function showExtenalLog(log_type)
 * 	->	This function is used to show in an external popup the log.
 * 
 * function fillUL(list, ul_id, icon)
 * 	->	Fill an UL element with items from a list items, should contains id and 
 * 		name properties.
 * 
 * function updateButtons(id, elm)
 * 	->	Function to update the buttons in the toolbar, after choose
 *  	a new protocol run.
 * 
 * function updateRow(id, elm, row)
 * 	->	Function to update the row in the protocol run list when an element is
 *		selected.
 * 
 * function updateTree(id, elm, row)
 * 	->	Function to update the node in the protocol run tree when an element is
 * 		selected.
 * 
 * function graphON(graph, icon_graph, list, icon_list)
 * 	->	Function to disable the list view and enable the tree graph view.
 * 
 * function graphOFF(graph, icon_graph, list, icon_list)
 * 	->	Function to disable the tree graph view and enable the list view.
 * 
 * function changeStatusGraph(status, graph, icon_graph, list, icon_list)
 * 	->	Function to switch between the graph/list view depending on the status.
 * 
 * function switchGraph() 
 * 	->	Main function called to change between graph tree/list views.
 * 
 * function updateGraphView(status)
 * 	->	Method to update the graph view in the internal settings of the project
 * 		passing like argument a boolean.
 * 		If the graph view was active when you closed the project last time will
 * 		be available directly, in the another case, will be inactive.
 * 
 * function editObject(objId)
 * 	->	Method to edit an object given his objId, calling an ajax to get the 
 * 		attributes from the object.
 * 
 * function deleteProtocolForm(protocolId)
 * 	->	Dialog web form based in messi.js to verify the option to delete a protocol. 
 * 
 * function deleteProtocol(elm)
 * 	->	Method to execute a delete for a protocol.
 * 
 * function copyProtocol(id)
 * 	-> Method to copy protocol.
 *
 * function stopProtocolForm(protocolId)
 * 	->	Dialog web form based in messi.js to verify the option to stop a protocol.
 * 
 * function stopProtocol(elm)
 * 	->	Method to stop the run for a protocol
 * 
 * function changeTreeView()
 * 	->	Method to update the protocol tree to run in the web left side.
 * 
 * function refreshRuns(mode)
 * 	->	Method to refresh the run list/graph checking changes. 
 * 
 * function checkStatusNode(id, status)
 * 	->	Function to check the status for a node (protocol) in the protocol graph.
 * 
 * 
 * function checkStatusRow(id, status, time)
 * 	-> Function to check the status for a node (protocol) in the protocol list.
 * 
 * 
 ******************************************************************************/

/** METHODS ***************************************************************** */

function isCtrlPress(event) {
	return event.ctrlKey;
}

function preloadToolbar(list) {
	/*
	 * Function to preload the nodes/rows selected
	 */
	if (list.length > 1) {
		refreshToolbarMultipleMark(list)
	} else {
		refreshToolbarSingleMark(list)
	}
}

function launchToolbarList(id, elm, multi) {
	/*
	 * Toolbar used in the project content template for list view
	 */

	if (multi) {
		enableMultipleMarkList(elm);
	} else {
		disableMultipleMarkList(id);
		launchToolbarProject(id, elm, "list");
	}
}

function launchToolbarTree(id, elm, multi) {
	/*
	 * Toolbar used in the project content template for graph view
	 */

	if (multi) {
		enableMultipleMarkGraph(elm);
	} else {
		disableMultipleMarkGraph(id);
		launchToolbarProject(id, elm, "graph");
	}
}

function launchToolbarProject(id, elm, type) {
	/*
	 * Method used to update an element in the template depending on view mode
	 * (list or graph)
	 */
	var row = $("div#toolbar");

	switch (type) {
	case "list":
		updateRow(id, elm, row);
		break;
	case "graph":
		updateTree(id, elm, row);
		break;
	}

	// Update the buttons functionalities into the toolbar
	// and update the content for the tabs
	refreshToolbarSingleMark(id)

}

function refreshToolbarSingleMark(id) {
	/*
	 * Update the buttons functionalities into the toolbar for single marks.
	 */
	updateButtons(id, "single");

	// Update the content for the tabs
	updateTabs(id);
}

function refreshToolbarMultipleMark(list_id) {
	/*
	 * Update the buttons functionalities into the toolbar for multiple marks.
	 */
	updateButtons(list_id, "multiple");
}

function getElmMarkedList() {
	/*
	 * Return a list of elements marked in the list view mode.
	 */
	var list_marked = [];

	$.each($("tr"), function() {
		var elm = $(this);
		var id = elm.attr("id");

		if (elm.hasClass("selected")) {
			list_marked.push(id);
		}
	});
	return list_marked;
}

function getElmMarkedGraph() {
	/*
	 * Return a list of elements marked in the graph view mode.
	 */
	var list_marked = [];

	$.each($(".window"), function() {
		var elm = $(this);
		var id = elm.attr("id").split("_").pop();

		if (elm.attr("selected") == "selected") {
			list_marked.push(id);
		}
	});
	return list_marked;
}

function transposeElmMarked(status) {
	/*
	 * Switch the elements marked between graph and list view.
	 */
	var list_marked = [];

	if (status == 'inactive') {
		// Transpose elements in the list marked to the graph

		$.each($("tr"), function() {
			var elm = $(this);
			var id = elm.attr("id");
			var elmToMark = $("#graph_" + id);
			var selected = elmToMark.attr("selected");

			if (elm.hasClass("selected")) {
				// console.log("adding "+ id)
				if (selected != "selected" || selected == undefined) {
					markSingleNodeGraph(elmToMark);
				}
				list_marked.push(id);
			} else if (!elm.hasClass("selected")
					&& elmToMark.attr("selected") == "selected") {
				dismarkSingleNodeGraph(elmToMark);
			}

		});
	} else if (status == 'active') {
		// Transpose elements in the graph marked to the list

		$.each($(".window"), function() {
			var elm = $(this);
			var id = elm.attr("id").split("_").pop();
			var elmToMark = $("tr#" + id);

			if (elm.attr("selected") == "selected") {
				// console.log("adding "+ id)
				if (!elmToMark.hasClass("selected")) {
					markSingleNodeList(elmToMark);
				}
				list_marked.push(id);
			} else if (elm.attr("selected") != "selected"
					&& elmToMark.hasClass("selected")) {
				dismarkSingleNodeList(elmToMark);
			}
		});
	}

	// console.log(list_marked)

	if (list_marked.length > 0) {
		// Update the runs selected in the DB
		refreshSelectedRuns(list_marked);
	}

}

function refreshSelectedRuns(list_marked) {
	/*
	 * Refresh the elements contained into the list.
	 */
	var URL = getSubDomainURL() + '/save_selection/?mark='
			+ listToString(list_marked)

	$.ajax({
		type : "GET",
		url : URL,
		async : false
	});
}

/** Graph Methods ********************************************** */

function markSingleNodeGraph(elm) {
	/*
	 * Method to mark a node selected in the graph.
	 */
	// Highlight the node
	elm.css("border", "2.5px solid Firebrick");
	elm.attr("selected", "selected");
}

function dismarkSingleNodeGraph(elm) {
	/*
	 * Method to dismark a node selected in the graph.
	 */

	// Clear the node
	elm.css("border", "");
	elm.removeAttr("selected")
}

function enableMultipleMarkGraph(elm) {
	/*
	 * Method to activate the multiple selection.
	 */
	if (elm.attr("selected") == "selected") {
		dismarkSingleNodeGraph(elm);
		$("div#graphActiv").removeAttr("data-option");
	} else {
		markSingleNodeGraph(elm);
		$("div#graphActiv").attr("data-option", elm.attr("id"));
	}
	var list_marked = getElmMarkedGraph()
	refreshToolbarMultipleMark(list_marked);
}

function disableMultipleMarkGraph(id) {
	/*
	 * Method to disactivate the multiple selection.
	 */
	$.each($(".window"), function() {
		var elm = $(this)
		if (elm.attr("id") != "graph_" + id
				&& elm.attr("selected") != undefined) {
			dismarkSingleNodeGraph(elm)
		}
	})
}

/** List Methods ********************************************** */

function markSingleNodeList(elm) {
	/*
	 * Method to mark a row selected in the list.
	 */

	// Highlight the node
	elm.addClass("selected");
}

function dismarkSingleNodeList(elm) {
	/*
	 * Method to dismark a row selected in the list.
	 */

	// Clear the node
	elm.removeClass("selected");
}

function enableMultipleMarkList(elm) {
	/*
	 * Method to activate the multiple selection.
	 */
	if (elm.hasClass("selected")) {
		dismarkSingleNodeList(elm);
	} else {
		markSingleNodeList(elm);
	}
	var list_marked = getElmMarkedList()
	refreshToolbarMultipleMark(list_marked);
}

function disableMultipleMarkList(id) {
	/*
	 * Method to disactivate the multiple selection.
	 */
	$.each($("tr"), function() {
		var elm = $(this)
		if (elm.hasClass("selected")) {
			dismarkSingleNodeList(elm);
		}
	})
}

/** *************************************************************************** */

function updateTabs(id) {
	/*
	 * Fill the content of the summary tab for a protocol run selected (Data /
	 * Summary / Methods / Status)
	 */
	var URL = getSubDomainURL() + '/protocol_info/?protocolId=' + id
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		success : function(json) {

			// DATA SUMMARY
			fillUL("input", json.inputs, "protocol_input", "fa-sign-in");
			fillUL("output", json.outputs, "protocol_output", "fa-sign-out");

			var summary = $("#protocol_summary");

			summary.empty();
			summary.append(json.summary);

			// METHODS
			$("#tab-methods").empty();
			$("#tab-methods").append(json.methods);

			// STATUS
			if (json.status == "running") {
				// Action Stop Button
				// $("span#analyzeTool").hide();
				$("span#buttonAnalyzeResult").hide();
				$("span#stopTool").show();
				$("a#stopTool").attr('href',
						'javascript:stopProtocolForm("' + id + '")');
			} else {
				// Action Analyze Result Button
				$("span#stopTool").hide();
				// $("span#analyzeTool").show();
				$("span#buttonAnalyzeResult").show();
				$("a#analyzeTool").attr('href',
						'javascript:launchViewer("' + id + '")');
				$("a#downloadTool").attr('href',
						'javascript:downloadOutput("' + id + '")');
			}

			// LOGS
			$("#tab-logs-output").empty();
			$("#tab-logs-output").append(json.logs_out);
			$("#tab-logs-error").empty();
			$("#tab-logs-error").append(json.logs_error);

			$("#tab-logs-scipion").empty();
			$("#tab-logs-scipion").append(json.logs_scipion);
		},
		error : function() {
			console.log("ERROR IN PROTOCOL_INFO REQUEST")
		}
	});
}

function showLog(log_type) {
	/*
	 * This function is used to show or hide differents logs about a protocol
	 * selected.
	 */

	switch (log_type) {
	case "output_log":

		$("div#tab-logs-output").css("display", "")
		$("a#output_log").attr("class", "elm-header-log_selected")
		html = $("div#tab-logs-output").html()

		$("div#tab-logs-error").css("display", "none")
		$("a#error_log").attr("class", "elm-header-log")

		$("div#tab-logs-scipion").css("display", "none")
		$("a#scipion_log").attr("class", "elm-header-log")

		break;

	case "error_log":

		$("div#tab-logs-output").css("display", "none")
		$("a#output_log").attr("class", "elm-header-log")

		$("div#tab-logs-error").css("display", "")
		$("a#error_log").attr("class", "elm-header-log_selected")
		html = $("div#tab-logs-error").html()

		$("div#tab-logs-scipion").css("display", "none")
		$("a#scipion_log").attr("class", "elm-header-log")

		break;

	case "scipion_log":

		$("div#tab-logs-output").css("display", "none")
		$("a#output_log").attr("class", "elm-header-log")

		$("div#tab-logs-error").css("display", "none")
		$("a#error_log").attr("class", "elm-header-log")

		$("div#tab-logs-scipion").css("display", "")
		$("a#scipion_log").attr("class", "elm-header-log_selected")
		html = $("div#tab-logs-scipion").html()

		break;
	}

	// Fill the button to show the button in an external window
	var log_func = "javascript:showExternalLog('" + log_type + "');"
	$("a#externalTool").attr("href", log_func)

}

function showExternalLog(log_id) {
	/*
	 * This function is used to show in an external popup the log.
	 */
	var html = "";

	switch (log_id) {

	case "output_log":
		html = $("div#tab-logs-output").html()
		break;
	case "error_log":
		html = $("div#tab-logs-error").html()
		break;
	case "scipion_log":
		html = $("div#tab-logs-scipion").html()
		break;
	}
	customPopupFileHTML(html, log_id, 1024, 768);

}

function fillUL(type, list, ulId, icon) {
	/*
	 * Fill an UL element with items from a list items should contains id and
	 * name properties
	 */
	var ul = $("#" + ulId);
	ul.empty();

	for (var i = 0; i < list.length; i++) {

		// Visualize Object
		var visualize_html = '<a href="javascript:launchViewer(' + list[i].id
				+ ');"><i class="fa ' + icon
				+ '" style="margin-right:10px;"></i>' + list[i].name

		var download_html = "";

		if (type == "input") {
			visualize_html += ' (from ' + list[i].nameId + ')</a>'
		} else if (type == "output") {
			visualize_html += '</a>'

			// Download File Object
			// download_html = '<a href="javascript:downloadOutput('+ list[i].id
			// + ');"> '+
			// '<i class="fa fa-save" style="margin-left:0px;">
			// Download</i></a>'

			// Update Tab Download Button
			// $("a#downloadTool").attr('href',
			// 'javascript:downloadOutput("'+list[i].id +'")');

		}

		// Edit Object
		var edit_html = '<a href="javascript:editObject(' + list[i].id
				+ ');"> '
				+ '<i class="fa fa-pencil" style="margin-left:0px;"></i></a>'

		ul.append('<li>' + visualize_html + edit_html + "&nbsp;&nbsp;&nbsp;"
				+ list[i].info + download_html + '</li>');

	}
}

function updateButtons(id, mode) {
	/*
	 * Function to update the buttons in the toolbar after choose a new protocol
	 * run.
	 */

	if (mode == undefined) {
		mode = "single";
	}

	switch (mode) {

		case "single":
	
			// Action Edit Button
			$("a#editTool").attr('href',
					'javascript:popup("/form/?protocolId=' + id + '")');
			$("span#editTool").show();
	
			// Action Copy Button
			$("a#copyTool").attr(
					'href',
					'javascript:popup("/form/?protocolId=' + id + '&action=copy'
							+ '")');
			$("span#copyTool").show();
	
			// Action Delete Button
			$("a#deleteTool").attr('href',
					'javascript:deleteProtocolForm("' + id + '","single")');
			$("span#deleteTool").show();
	
			// Action Browse Button
			var aux = "javascript:alert('Not implemented yet')";
			$("a#browseTool").attr('href', aux);
			$("span#browseTool").show();
	
			break;
	
		case "multiple":
	
			// Action Edit Button
			$("a#editTool").attr('href', '#');
			$("span#editTool").hide();
	
			// Action Copy Button
			$("a#copyTool").attr('href', 'javascript:copyProtocol("' + id + '")');
			$("span#copyTool").show();
	
			// Action Delete Button
			$("a#deleteTool").attr('href',
					'javascript:deleteProtocolForm("' + id + '","multiple")');
			$("span#deleteTool").show();
	
			// Action Browse Button
			$("a#browseTool").attr('href', '#');
			$("span#browseTool").hide();
	
			break;

	}
	// Show toolbar
	$("div#toolbar").show();
}

function updateRow(id, elm, row) {
	/*
	 * Function to update the row in the protocol run list when an element is
	 * selected.
	 */
	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
		dismarkSingleNodeList(rowOld);
	}
	markSingleNodeList(elm)

	// add id value into the toolbar
	row.attr('value', id);
}

function updateTree(id, elm, row) {
	/*
	 * Function to update the node in the protocol run tree when an element is
	 * selected.
	 */
	var oldSelect = $("div#graphActiv").attr("data-option");
	var selected = "graph_" + id;

	if (oldSelect != selected) {
		if (oldSelect != "") {
			var aux = "div#" + oldSelect + ".window";
			dismarkSingleNodeGraph($(aux))
		}
		row.attr('value', id);
		$("div#graphActiv").attr("data-option", selected);
		markSingleNodeGraph(elm)
	}
}

function switchView(state, elm, tool){
	if (state){
		elm.attr("data-mode", "active");
		elm.show();
		tool.hide();
	} else{
		elm.attr("data-mode", "inactive");
		elm.hide();
		tool.show();
	}
}

function changeStatusGraph(mode, param) {
	/*
	 * Function to switch between the graph/list/small view depending on the status.
	 */

	graph = param['graph'];
	graphTool = param['graphTool'];
	list = param['list'];
	listTool = param['listTool'];
	smallGraph = param['smallGraph'];
	smallGraphTool = param['smallGraphTool'];

	if (mode == 1) {
		switchView(true, graph, graphTool)
		switchView(false, list, listTool)
		switchView(false, smallGraph, smallGraphTool)
		updateGraphView("1");
	} else if (mode == 0) {
		switchView(false, graph, graphTool)
		switchView(true, list, listTool)
		switchView(false, smallGraph, smallGraphTool)
		updateGraphView("0");
	} else if (mode == 2) {
		switchView(false, graph, graphTool)
		switchView(false, list, listTool)
		switchView(true, smallGraph, smallGraphTool)
		updateGraphView("2");
	}
}

function switchMode(mode) {
	/*
	 * Main function called to change between graph tree/list/small views.
	 */

	// element marked obtained from value in the toolbar
	var id = $("div#toolbar").attr("value");
	
	// get row elements
	var graph = $("div#graphActiv");
	var graphTool = $("span#treeTool");
	var list = $("div#runTable");
	var listTool = $("span#listTool");
	var smallGraph = $("div#graphSmallActiv");
	var smallGraphTool = $("span#treeSmallTool");

	// create dictionary of elements
	var param = {};
	param['graph'] = graph;
	param['graphTool'] = graphTool; 
	param['list'] = list;
	param['listTool'] = listTool;
	param['smallGraph'] = smallGraph;
	param['smallGraphTool'] = smallGraphTool ;
	
	changeStatusGraph(mode, param);
	
	if (graph.attr("data-time") == 'first' && mode == 1) {
		// Graph will be painted once
		callPaintGraph(graph, "normal");
		graph.attr("data-time", "not");
	}
	
	if (smallGraph.attr("data-time") == 'first' && mode == 2) {
		// Small Graph will be painted once
		callPaintGraph(smallGraph, "small");
		smallGraph.attr("data-time", "not");
	}

	// Keep the consistency about the selected elements between
	// the list and graph views.
	
	if (mode == 1){
		transposeElmMarked("active");
	}
	else {
		transposeElmMarked("inactive");
	}

}

function updateGraphView(status) {
	/*
	 * Method to update the graph view in the internal settings of the project
	 * passing like argument a boolean. If the graph view was active when you
	 * closed the project last time will be available directly, in the another
	 * case, will be inactive.
	 */
	var URL = getSubDomainURL() + "/update_graph_view/?status=" + status
	$.ajax({
		type : "GET",
		url : URL,
		async : false
	});
}

function editObject(objId) {
	var URL = getSubDomainURL() + '/get_attributes/?objId=' + objId
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "text",
		success : function(text) {
			var res = text.split("_-_")
			label = res[0]
			comment = res[1]
			editObjParam(objId, "Label", label, "Comment", comment,
					"Describe your run here...", "object")
		}
	});
}

function downloadOutput(objId) {
	var URL = getSubDomainURL() + '/download_output/?objId=' + objId
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "text",
		success : function(text) {
			var URL = getSubDomainURL() + '/get_file/?path=' + text
					+ '&filename=output.zip'
			window.location.href = URL;
		}
	});
}

function deleteProtocolForm(id, mode) {
	/*
	 * Dialog web form based in messi.js to verify the option to delete.
	 */

	var msg = "</td><td class='content' value='" + id

	switch (mode) {
	case "single":
		msg += "'><strong>ALL DATA</strong> related to this <strong>protocol run</strong>"
		break;
	case "multiple":
		msg += "'><strong>ALL DATA</strong> related to the <strong>selected protocols run</strong>"
		break;
	}
	msg += " will be <strong>DELETED</strong>. Do you really want to continue?</td></tr></table>";

	warningPopup('Confirm DELETE', msg, 'deleteProtocol')
}

function deleteProtocol(elm) {
	/*
	 * Method to execute a delete for a protocol
	 */
	var id = elm.attr('value');
	var URL = getSubDomainURL() + "/delete_protocol/?id=" + id
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		async : false,
		success : function(json) {
			if (json.errors != undefined) {
				// Show errors in the validation
				errorPopup('Errors found', json.errors);
			} else if (json.success != undefined) {
				// launchMessiSimple("Successful", messiInfo(json.success));
				// window.location.reload()
			}
		},
		error : function() {
			alert("error")
		}
	});

}

function copyProtocol(id) {
	/*
	 * Method to copy protocol.
	 */
	var URL = getSubDomainURL() + "/copy_protocol/?id=" + id
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		async : false,
		success : function(json) {
			if (json.errors != undefined) {
				// Show errors in the validation
				errorPopup('Errors found', json.errors);
			} else if (json.success != undefined) {
				// launchMessiSimple("Successful", messiInfo(json.success));
				// window.location.reload()
			}
		},
		error : function() {
			alert("error")
		}
	});
}

function stopProtocolForm(protocolId) {
	/*
	 * Dialog web form based in messi.js to verify the option to stop a
	 * protocol.
	 */

	var msg = "<td class='content' value='"
			+ protocolId
			+ "'>This <strong>protocol run</strong>"
			+ " will be <strong>STOPPED</strong>. Do you really want to continue?</td>";

	warningPopup('Confirm STOP', msg, 'stopProtocol');
}

function stopProtocol(elm) {
	/*
	 * Method to stop the run for a protocol
	 */
	var protId = elm.attr('value');
	var URL = getSubDomainURL() + "/stop_protocol/?protocolId=" + protId
	$.ajax({
		type : "GET",
		url : URL
	});
}

function changeTreeView() {
	/*
	 * Method to update the protocol tree to run in the web left side.
	 */
	protIndex = $('#viewsTree').val();
	var URL = getSubDomainURL() + '/update_prot_tree/?index=' + protIndex

	$.ajax({
		type : "GET",
		url : URL,
		dataType : "text",
		success : function() {
			var URL2 = getSubDomainURL() + '/tree_prot_view/'
			$.ajax({
				type : "GET",
				url : URL2,
				dataType : "html",
				async : false,
				success : function(data) {
					$('div.protFieldsetTree').html(data);
				}
			});
		}
	});
}

// ** REFRESHING FUNCTIONALITIES
// *****************************************************/

function refreshRuns(mode) {
	/*
	 * Method to update the run list/graph
	 */

	var URL = getSubDomainURL() + '/run_table_graph/'

	$(function() {
		$.ajax({
			async : true,
			url : URL,
			// datatype: "text",
			success : function(data) {

//				console.log("data: " + data)

				if (typeof data == 'string' || data instanceof String) {

					if (data == 'stop') {
						window.clearTimeout(updatetimer);
						// stop the script
					} else if (data == 'ok') {
						// no changes
					} else {
						$('div#runsInfo').html(data);
						// refresh the data keeping the element marked
					}
				} else {
					for (var x = 0; x < data.length; x++) {
						var id = data[x][0];
						var status = data[x][2];
						var time = data[x][3];

						checkStatusNode(id, status)
						checkStatusRow(id, status, time)

//						console.log(status)

						if (status == "finished") {
//							console.log("protocol finished!")
							// updateTabs(id)
						}
					}
				}
			},
			error : function() {
				console.log("Error in the refresh")
			}
		});
	});

	if (mode) {
		var updatetimer = setTimeout(function() {
			refreshRuns(1);
		}, 3000);
	}
}

function checkStatusNode(id, status) {
	node = $("#graph_" + id);

	if (status.indexOf('running') == 0) {
		node.css("background-color", "#FCCE62");
	}
	if (status.indexOf('failed') == 0) {
		node.css("background-color", "#F5CCCB");
	}
	if (status.indexOf('finished') == 0) {
		node.css("background-color", "#D2F5CB");
	}
	node.find("#nodeStatus").html(status);
}

function checkStatusRow(id, status, time) {
	row = $("tr#" + id);
	row.find(".status").html(status);
	row.find(".time").html(time);
}

function openSearchProtocolPopup(search)
{
	$.ajax({
		type : "GET",
		url : getSubDomainURL() + "/search_protocol/",
		success : function(html) {
			new Messi(html, 
			{
				title : 'Search protocol',
				modal : true,
				buttons : [ {
					id : 1,
					label : 'Ok',
					val : 'C'
				}]
			});
		}
	});
	    
		
	
}




