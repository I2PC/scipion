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
 *  e-mail address 'jmdelarosa@cnb.csic.es'
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
 * function launchToolbarList(id, elm)
 * 	->	Toolbar used in the project content template for the run list view.
 * 
 * function launchToolbarList(id, elm)
 * 	->	Toolbar used in the project content template for the runs tree view.
 * 
 * function checkRunStatus(id)
 * 	->	Function to check a protocol run, depend on the status two button will be
 * 		switching in the toolbar (Stop / Analyze Results).
 * 
 * function fillTabsSummary(id)
 * 	->	Fill the content of the summary tab for a protocol run selected 
 * 		(Data and Summary)
 * 
 * function fillUL(list, ul_id, icon)
 * 	->	Fill an UL element with items from a list items, should contains id and 
 * 		name properties.
 * 
 * function launchViewer(id)
 * 	->	Launch the viewers to analyze the results for a protocol run.
 * 
 * function updateButtons(id, elm)
 * 	->	Function to update the buttons in the toolbar and the tabs, after choose
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
 * function markElmGraph(node_id, graph)
 * 	->	Function used to mark the same protocol run when one is selected in the
 * 		protocol run list and the view is changed to protocol run graph.
 * 
 * function markElmList(row_id, graph)
 * 	->	Function used to mark the same protocol run when one is selected in the
 * 		protocol run graph and the view is changed to protocol run list.
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
 * function deleteProtocolForm(protocolId)
 * 	->	Dialog web form based in messi.js to verify the option to delete a protocol. 
 * 
 * function deleteProtocol(elm)
 * 	->	Method to execute a delete for a protocol.
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
 ******************************************************************************/

 /** METHODS ******************************************************************/

function launchToolbarList(id, elm) {
	/*
	 * Toolbar used in the project content template for list view
	 */
	var row = $("div#toolbar");
	updateRow(id, elm, row);
	updateButtons(id, elm);
	row.show(); // Show toolbar
}

function launchToolbarTree(id, elm) {
	/*
	 * Toolbar used in the project content template for list view
	 */
	var row = $("div#toolbar");
	updateTree(id, elm, row);
	updateButtons(id, elm);
	row.show(); // Show toolbar
}

function checkRunStatus(id) {
	/*
	 * Function to check a protocol run, depend on the status two button will be
	 * switching in the toolbar (Stop / Analyze Results).
	 */
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
				'javascript:stopProtocolForm("' + id + '")');
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

function fillTabsSummary(id) {
	/*
	 * Fill the content of the summary tab for a protocol run selected 
	 * (Data and Summary)
	 */
	$.ajax({
		type : "GET",
		url : '/protocol_io/?protocolId=' + id,
		dataType : "json",
		success : function(json) {
			fillUL(json.inputs, "protocol_input", "fa-sign-in");
			fillUL(json.outputs, "protocol_output", "fa-sign-out");
		}
	});

	$.ajax({
		type : "GET",
		url : '/protocol_summary/?protocolId='+ id,
		dataType : "json",
		success : function(json) {
			$("#tab-summary").empty();
			$("#tab-summary").append(json);
//			for ( var i = 0; i < json.length; i++) {
//				$("#tab-summary").append('<p>' + json[i] + '</p>');
//			}
		}
	});
}

function fillUL(list, ulId, icon) {
	/*
	 * Fill an UL element with items from a list items should contains id and name
	 * properties
	 */
	var ul = $("#" + ulId);
	ul.empty();
	for ( var i = 0; i < list.length; i++) {
//		ul.append('<li><a href="/visualize_object/?objectId=' + list[i].id
//				+ '"target="_blank"><img src="../../../../resources/' + icon + '" /> '
//				+ list[i].name + '</a></li>');
		ul.append('<li><a href="/visualize_object/?objectId=' + list[i].id
				+ '"target="_blank"><i class="fa ' + icon + '" style="margin-right:10px;"></i>'
				+ list[i].name + '</a></li>');
	}
}

function launchViewer(id){
	/*
	 * Launch the viewers to analyze the results of the protocol run
	 */
	$.ajax({
		type : "GET",
		// Execute the viewer 
		url : "/launch_viewer/?protocolId=" + id,
		dataType : "json",
		success : function(json) {
			popUpJSON(json);
		}
	});	
}

function updateButtons(id, elm){
	/*
	 * Function to update the buttons in the toolbar and the tabs, after choose a new protocol run.
	 */
	// Action Edit Button
	$("a#editTool").attr('href',
	'javascript:popup("/form/?protocolId=' + id + '")');
	
	// Action Copy Button
	$("a#copyTool").attr('href',
	'javascript:popup("/form/?protocolId=' + id + '&action=copy' + '")');

	// Action Delete Button
	$("a#deleteTool").attr('href',
			'javascript:deleteProtocolForm("' + id + '")');

	// Action Browse Button
	var aux = "javascript:alert('Not implemented yet')";
	$("a#browseTool").attr('href', aux);
	
	checkRunStatus(id);
	fillTabsSummary(id);
}

function updateRow(id, elm, row){	
	/*
	 * Function to update the row in the protocol run list when an element is
	 * selected.
	 */
	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.removeClass('selected')
	}
	elm.addClass('selected')

	// add id value into the toolbar
	row.attr('value', id);
}

function updateTree(id, elm, row){
	/*
	 * Function to update the node in the protocol run tree when an element is
	 * selected.
	 */
	var oldSelect = $("div#graphActiv").attr("data-option");
	var selected = "graph_" + id;

	if (oldSelect != selected) {
		if (oldSelect != "") {
			var aux = "div#" + oldSelect + ".window";
			$(aux).css("border", "");
		}
		row.attr('value', id);
		$("div#graphActiv").attr("data-option", selected);
		elm.css("border", "2.5px solid Firebrick");
	}
}

function graphON(graph, icon_graph, list, icon_list){
	/*
	 * Function to disable the list view and enable the tree graph view.
	 */
	
	// Graph ON
	graph.attr("data-mode", "active");
	graph.attr("style", "");
	icon_graph.hide();
	
	// Table OFF
	list.attr("data-mode", "inactive");
	list.attr("style", "display:none;");
	icon_list.show();

	// Update Graph View
	updateGraphView("True");
}

function graphOFF(graph, icon_graph, list, icon_list){
	/*
	 * Function to disable the tree graph view and enable the list view.
	 */
	
	// Table ON	
	list.attr("data-mode", "active");
	list.attr("style", "");
	icon_list.hide();
	
	// Graph OFF
	graph.attr("data-mode", "inactive");
	graph.attr("style", "display:none;");
	icon_graph.show();
	
	// Update Graph View
	updateGraphView("False")
}

function changeStatusGraph(status, graph, graphTool, list, listTool){
	/*
	 * Function to switch between the graph/list view depending on the status.
	 */
	if (status == 'inactive') {
		// Graph ON & Table OFF
		graphON(graph, graphTool, list, listTool);
	} else if (status == 'active') {
		// Table ON	& Graph OFF
		graphOFF(graph, graphTool, list, listTool);
	}
}

function markElmGraph(node_id, graph){
	/*
	 * Function used to mark the same protocol run when one is selected in the
	 * protocol run list and the view is changed to protocol run tree.
	 */
	var s = "graph_" + node_id;

	if (s != "" || s != undefined) {
		var nodeClear = graph.attr("data-option");
		
		if (nodeClear.length>0 && nodeClear != undefined) {
			// Clear the node
			var elmClear = $("div#" + nodeClear);
			elmClear.css("border", "");
		} 
		// setElement in graph
		graph.attr("data-option", s);
	
		// Highlight the node
		var elm = $("div#" + s);
		elm.css("border", "2.5px solid Firebrick");
	}
}
	
function markElmList(row_id, graph){
	/*
	 * Function used to mark the same protocol run when one is selected in the
	 * protocol run tree and the view is changed to protocol run list.
	 */
	var rowClear = $("tr.selected").attr("id");
	if (rowClear != "") {
		if (rowClear != row_id) {
			// Clear the row selected
			var elmClear = $("tr.selected");
			elmClear.attr("style", "");
			elmClear.attr("class", "runtr");

			// setElement in table
			var elm = $("tr#" + row_id + ".runtr");
			var projName = graph.attr("data-project");
			launchToolbarList(row_id, elm);
		}
	}
}

function switchGraph() {
	/*
	 * Main function called to change between graph tree/list views.
	 */
	
	// graph status (active or inactive) 
	var status = $("div#graphActiv").attr("data-mode");

	// element marked obtained from value in the toolbar
	var id = $("div#toolbar").attr("value");
	
	//get row elements 
	var graph = $("div#graphActiv");
	var graphTool = $("span#treeTool");
	var list = $("div#runTable");
	var listTool = $("span#listTool");
	
	changeStatusGraph(status, graph, graphTool, list, listTool)
		
	// Graph will be painted once
	if (graph.attr("data-time") == 'first') {
		callPaintGraph();
		graph.attr("data-time", "not");
	} 
	markElmGraph(id, graph);
	markElmList(id, graph);
}

function updateGraphView(status) {
	/*
	 * Method to update the graph view in the internal settings of the project
	 * passing like argument a boolean.
	 * If the graph view was active when you closed the project last time will
	 * be available directly, in the another case, will be inactive.
	 */
	$.ajax({
		type : "GET", 
		url : "/update_graph_view/?status=" + status
	});
}

function deleteProtocolForm(protocolId) {
	/*
	 * Dialog web form based in messi.js to verify the option to delete.
	 */
	var msg = "</td><td class='content' value='"
			+ protocolId
			+ "'><strong>ALL DATA</strong> related to this <strong>protocol run</strong>"
			+ " will be <strong>DELETED</strong>. Do you really want to continue?</td></tr></table>";
	
	msg = messiWarning(msg);

	new Messi(msg, {
		title : 'Confirm DELETE',
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : 'deleteProtocol'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'fa-ban'
		} ]
//		callback : function(val) {
//			if (val == 'Y') {
//				window.location.href = "/project_content/?projectName="
//						+ projName;
//			}
//		}
	});
}

function deleteProtocol(elm) {
	/*
	 * Method to execute a delete for a protocol
	 */
	var protId = elm.attr('value');
	
	$.ajax({
		type : "GET",
		url : "/delete_protocol/?protocolId=" + protId,
		dataType : "json",
		success : function(json) {
			if(json.errors != undefined){
				// Show errors in the validation
				showErrorValidation(json.errors);
			} else if(json.success!= undefined){
//				launchMessiSimple("Successful", messiInfo(json.success));
//				window.location.reload()
			}
		},
		error: function(){
			alert("error")
		}
	});
}

function stopProtocolForm(protocolId) {
	/*
	 * Dialog web form based in messi.js to verify the option to stop a protocol.
	 */
		
	var msg = "<td class='content' value='"
			+ protocolId
			+ "'>This <strong>protocol run</strong>"
			+ " will be <strong>STOPPED</strong>. Do you really want to continue?</td>";
			
	msg = messiWarning(msg);

	new Messi(msg, {
		title : 'Confirm STOP',
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : 'stopProtocol'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'fa-ban'
		} ]
//		callback : function(val) {
//			if (val == 'Y') {
//				window.location.href = "/project_content/?projectName="
//						+ projName;
//			}
//		}
	});
}

function stopProtocol(elm) {
/*
 * Method to stop the run for a protocol
 */
	var protId = elm.attr('value');
	
	$.ajax({
		type : "GET",
		url : "/stop_protocol/?protocolId=" + protId
	});
}

function changeTreeView(){
	/*
	 * Method to update the protocol tree to run in the web left side.
	 */
	protIndex = $('#viewsTree').val();
	
	$.ajax({
		type : "GET",
		url : '/update_prot_tree/?index='+ protIndex,
		dataType:"text",
		success : function() {
			$.ajax({
				url: '/tree_prot_view/',
				success: function(data) {
					$('div.protFieldsetTree').html(data);
				}
			});
		}
	});
}

function refreshRuns(mode){
	/*
	 * Method to update the run list/graph
	 */
	
	$(function() {
		$.ajax({
			url : '/run_table_graph/',
			success : function(data) {
				if (data=='stop'){
					window.clearTimeout(updatetimer);
					// stop the script
				}
				else if(data == 'ok'){
					// no changes
				}
				else {
					$('div#runsInfo').html(data);
					// refresh the data keeping the element marked
				}
			}
		});
  	});
	
	if(mode){
		var updatetimer = setTimeout(function(){ 
			refreshRuns(1);
	  	}, 3000);
	}
}



