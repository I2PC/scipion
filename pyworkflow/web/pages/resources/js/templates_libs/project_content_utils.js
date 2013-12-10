 /**************************************************************************
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
 **************************************************************************/



/**
 * Methods used in the project_content template. 
 * Toolbar functions + Manage tabs
 * 
 * launchToolbarList(projName, id, elm); 
 * launchToolbarTree(projName, id, elm);
 * checkRunStatus(projName, id);
 * fillTabsSummary(id);
 * fillUL(list, ulId, icon);
 * launchViewer(id);
 * updateButtons(projName, id, elm)
 * updateTree(id, elm);
 * updateRow(id, elm, row);
 * switchGraph();
 * deleteProtocolForm(projName, protocolId);
 * deleteProtocol(elm);
 * stopProtocolForm(projName, protocolId);
 * stopProtocol(elm);
 * 
 **/

/*
 * Toolbar used in the project content template for list view
 */
function launchToolbarList(id, elm) {
	var row = $("div#toolbar");
	updateRow(id, elm, row);
	updateButtons(id, elm);
	row.show(); // Show toolbar
}

/*
 * Toolbar used in the project content template for list view
 */
function launchToolbarTree(id, elm) {
	var row = $("div#toolbar");
	updateTree(id, elm, row);
	updateButtons(id, elm);
	row.show(); // Show toolbar
}

function checkRunStatus(id) {
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

/*
 * Fill the content of the summary tab
 */
function fillTabsSummary(id) {
	$.ajax({
		type : "GET",
		url : '/protocol_io/?protocolId=' + id,
		dataType : "json",
		success : function(json) {
//			fillUL(json.inputs, "protocol_input", "db_input.gif");
			fillUL(json.inputs, "protocol_input", "fa-sign-in");
//			fillUL(json.outputs, "protocol_output", "db_output.gif");
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

/*
 * Fill an UL element with items from a list items should contains id and name
 * properties
 */
function fillUL(list, ulId, icon) {
	var ul = $("#" + ulId);
	ul.empty();
	for ( var i = 0; i < list.length; i++) {
//		ul.append('<li><a href="/visualize_object/?objectId=' + list[i].id
//				+ '"target="_blank"><img src="../../../../resources/' + icon + '" /> '
//				+ list[i].name + '</a></li>');
		ul.append('<li><a href="/visualize_object/?objectId=' + list[i].id
				+ '"target="_blank"><i class="fa ' + icon + '"></i>'
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
		url : "/launch_viewer/?protocolId=" + id,
		dataType : "json",
		success : function(json) {
			popUpJSON(json);
		}
	});	
}

function updateButtons(id, elm){
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

function updateTree(id, elm, row){
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
	
function updateRow(id, elm, row){	
	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.removeClass('selected')
	}
	elm.addClass('selected')

	// add id value into the toolbar
	row.attr('value', id);
}

function graphON(graph, graphTool, list, listTool){
	// Graph ON
	graph.attr("data-mode", "active");
	graph.attr("style", "");
	graphTool.hide();
	
	// Table OFF
	list.attr("data-mode", "inactive");
	list.attr("style", "display:none;");
	listTool.show();

	// Update Graph View
	updateGraphView("True");
}

function graphOFF(graph, graphTool, list, listTool){
	// Table ON	
	list.attr("data-mode", "active");
	list.attr("style", "");
	listTool.hide();
	// Graph OFF
	graph.attr("data-mode", "inactive");
	graph.attr("style", "display:none;");
	graphTool.show();
	
	// Update Graph View
	updateGraphView("False")
}

function changeStatusGraph(status, graph, graphTool, list, listTool){
	if (status == 'inactive') {
		// Graph ON & Table OFF
		graphON(graph, graphTool, list, listTool);
	} else if (status == 'active') {
		// Table ON	& Graph OFF
		graphOFF(graph, graphTool, list, listTool);
	}
}

function markElmGraph(id, graph){
	var s = "graph_" + id;

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
	
function markElmList(id, graph){
	
	var rowClear = $("tr.selected").attr("id");
	if (rowClear != "") {
		if (rowClear != id) {
			// Clear the row selected
			var elmClear = $("tr.selected");
			elmClear.attr("style", "");
			elmClear.attr("class", "runtr");

			// setElement in table
			var elm = $("tr#" + id + ".runtr");
			var projName = graph.attr("data-project");
			launchToolbarList(id, elm);
		}
	}
}


function switchGraph() {
	// graph status (active or inactive) 
	var status = $("div#graphActiv").attr("data-mode");

	// element marked obtained from value in the toolbar
	var id = $("div#toolbar").attr("value");
	
	//get row elements 
	var graph = $("div#graphActiv");
	var graphTool = $("div#treeTool");
	var list = $("div#runTable");
	var listTool = $("div#listTool");
	
	changeStatusGraph(status, graph, graphTool, list, listTool)
		
	// Graph will be painted once
	if (graph.attr("data-time") == 'first') {
		callPaintGraph();
		graph.attr("data-time", "not");
	} 
	
	markElmGraph(id, graph)
	markElmList(id, graph)
	
}


function updateGraphView(status) {
	$.ajax({
		type : "GET", 
		url : "/update_graph_view/?status=" + status
	});
}

/*
 * Dialog form to verify the right option to delete
 */
function deleteProtocolForm(protocolId) {

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
			btnClass : 'btn-select',
			btnFunc : 'deleteProtocol'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'btn-cancel'
		} ]
//		callback : function(val) {
//			if (val == 'Y') {
//				window.location.href = "/project_content/?projectName="
//						+ projName;
//			}
//		}
	});
}

/*
 * Method to execute a delete by a protocol
 */
function deleteProtocol(elm) {
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

/*
 * Dialog form to verify the right option to stop a protocol
 */
function stopProtocolForm(protocolId) {
		
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
			btnClass : 'btn-select',
			btnFunc : 'stopProtocol'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'btn-cancel'
		} ]
//		callback : function(val) {
//			if (val == 'Y') {
//				window.location.href = "/project_content/?projectName="
//						+ projName;
//			}
//		}
	});
}

/*
 * Method to stop the run for a protocol
 */
function stopProtocol(elm) {
	var protId = elm.attr('value');
	
	$.ajax({
		type : "GET",
		url : "/stop_protocol/?protocolId=" + protId
	});
}

/*
 * Method to update the protocol tree
 */
function changeTreeView(){
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

/*
 * Method to update the run list/graph
 */
function refreshRuns(){
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
	
	var updatetimer = setTimeout(function(){ 
		refreshRuns();
  	}, 3000);
}

