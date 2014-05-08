 /*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
 * 			   Adrian Quintana (aquintana@cnb.csic.es)
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
 * function enableMultipleMarkGraph(elm)
 * 	->	Used to highlight a node in the protocol graph.
 * 		This method is used to multiples marked nodes.
 * 
 * function disableMultipleMarkGraph(elm)
 * 	->	Used to remove the highlight applied to a node in the protocol graph.
 * 		This method is used to multiples marked nodes.
 * 
 * function launchToolbarTree(id, elm)
 * 	->	Toolbar used in the project content template for the runs tree view.
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
 * function markElmGraphSingle(node_id, graph)
 * 	->	Function used to mark the same protocol run when one is selected in the
 * 		protocol run list and the view is changed to protocol run graph.
 * 
 * function markElmListSingle(row_id, graph)
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

 /** METHODS ******************************************************************/

/** EVENTS WITH TRIGGERS **/
var event_graph = jQuery("div#graphActiv").trigger(jQuery.Event("click"));
var event_list = jQuery("div#runTable").trigger(jQuery.Event("click"));
var event = event_graph || event_list ; 


function launchToolbarList(id, elm) {
	/*
	 * Toolbar used in the project content template for list view
	 */
	if (event.ctrlKey){
		enableMultipleMarkList(elm);
	} else {
		disableMultipleMarkList(id);
		launchToolbarProject(id, elm, "list");
	}
}


function launchToolbarTree(id, elm) {
	/*
	 * Toolbar used in the project content template for graph view
	 */
	if (event.ctrlKey){
		enableMultipleMarkGraph(elm);
	} else {
		disableMultipleMarkGraph(id);
		launchToolbarProject(id, elm, "graph");
	}
}

	
function launchToolbarProject(id, elm, type){
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
	refreshToolbarSingleMark(id, elm)
	
	// Update the content for the tabs
	updateTabs(id);
}

function refreshToolbarSingleMark(id, elm){
	// Update the buttons functionalities into the toolbar
	updateButtons(id, "single");
}

function refreshToolbarMultipleMark(list_id){
	// Update the buttons functionalities into the toolbar
	updateButtons(list_id, "multiple");
}

function getElmMarkedList(){
	var list_marked = [];

	$.each($("tr"), function(){
		var elm = $(this);
		var id = elm.attr("id");
		
		if(elm.hasClass("selected")){
			list_marked.push(id);
		}
	});
	return list_marked;
}

function getElmMarkedGraph(){
	var list_marked = [];

	$.each($(".window"), function(){
		var elm = $(this);
		var id = elm.attr("id").split("_").pop();
		
		if(elm.attr("selected") == "selected"){
			list_marked.push(id);
		} 
	});
	return list_marked;
}


function transposeElmMarked(status){
	var list_marked = [];

	if (status == 'inactive') {
		// Transpose elements in the list marked to the graph
		
		$.each($("tr"), function(){
			var elm = $(this);
			var id = elm.attr("id");
			var elmToMark = $("#graph_"+id);
			var selected = elmToMark.attr("selected"); 
			
			if(elm.hasClass("selected") && (selected != "selected" || selected == undefined)){
//				console.log("adding "+ id)
				markSingleNodeGraph(elmToMark);
				list_marked.push(id);
			} else if(!elm.hasClass("selected") && elmToMark.attr("selected")=="selected"){
				dismarkSingleNodeGraph(elmToMark);
			} 
			
		});
	}
	else if (status == 'active') {
		// Transpose elements in the graph marked to the list
		
		$.each($(".window"), function(){
			var elm = $(this);
			var id = elm.attr("id").split("_").pop();
			var elmToMark = $("tr#"+id);
			
			if(elm.attr("selected") == "selected" && !elmToMark.hasClass("selected")){
//				console.log("adding "+ id)
				markSingleNodeList(elmToMark);
				list_marked.push(id);
			} else if (elm.attr("selected") != "selected" && elmToMark.hasClass("selected")){
				dismarkSingleNodeList(elmToMark);
			}
		});
	}
	
	if (list_marked.length > 0){
		// Update the runs selected in the DB
		refreshSelectedRuns(list_marked);
	}
	
}

function refreshSelectedRuns(list_marked){
	$.ajax({
		type : "GET", 
		url : '/save_selection/?mark=' + listToString(list_marked),
		async: false
	});
}

	
/** Graph Methods ***********************************************/

function markSingleNodeGraph(elm){
	// Highlight the node
	elm.css("border", "2.5px solid Firebrick");
	elm.attr("selected", "selected");
}

function dismarkSingleNodeGraph(elm){
	// Clear the node
	elm.css("border", "");
	elm.removeAttr("selected")
}

function enableMultipleMarkGraph(elm){
	if (elm.attr("selected") == "selected"){
		dismarkSingleNodeGraph(elm);
		$("div#graphActiv").removeAttr("data-option");
	} else {
		markSingleNodeGraph(elm);
		$("div#graphActiv").attr("data-option", elm.attr("id"));
	}
	var list_marked = getElmMarkedGraph()
	refreshToolbarMultipleMark(list_marked);
} 	

function disableMultipleMarkGraph(id){
	$.each($(".window"), function(){
		var elm = $(this)
		if (elm.attr("id") != "graph_"+id && elm.attr("selected") != undefined){
			dismarkSingleNodeGraph(elm)
		}
	}) 
}

/** List Methods ***********************************************/

function markSingleNodeList(elm){
	// Highlight the node
	elm.addClass("selected");
}

function dismarkSingleNodeList(elm){
	// Clear the node
	elm.removeClass("selected");
}

function enableMultipleMarkList(elm){
	if (elm.hasClass("selected")){
		dismarkSingleNodeList(elm);
	} else {
		markSingleNodeList(elm);
	}
	var list_marked = getElmMarkedList()
	refreshToolbarMultipleMark(list_marked);
}

function disableMultipleMarkList(id){
	$.each($("tr"), function(){
		var elm = $(this)
		if(elm.hasClass("selected")){
			dismarkSingleNodeList(elm);
		}
	}) 
}
	
/******************************************************************************/


function updateTabs(id) {
	/*
	 * Fill the content of the summary tab for a protocol run selected 
	 * (Data / Summary / Methods / Status)
	 */
	$.ajax({
		type : "GET",
		url : '/protocol_info/?protocolId=' + id,
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
			if(json.status=="running"){
				// Action Stop Button
//				$("span#analyzeTool").hide();
				$("span#buttonAnalyzeResult").hide();
				$("span#stopTool").show();
				$("a#stopTool").attr('href',
				'javascript:stopProtocolForm("' + id + '")');
			} else {
				// Action Analyze Result Button
				$("span#stopTool").hide();
//				$("span#analyzeTool").show();
				$("span#buttonAnalyzeResult").show();
				$("a#analyzeTool").attr('href', 'javascript:launchViewer("'+id +'")');
			}
			
			//LOGS
			$("#tab-logs-output").empty();
			$("#tab-logs-output").append(json.logs_out);
			
			$("#tab-logs-error").empty();
			$("#tab-logs-error").append(json.logs_error);
			
			$("#tab-logs-scipion").empty();
			$("#tab-logs-scipion").append(json.logs_scipion);
		}
	});
}

function showLog(log_type){
	
	switch (log_type) {
		case "output_log":
			
			$("#tab-logs-output").css("display","")
			$("#output_log").attr("class", "elm-header-log_selected")
			html = $("#tab-logs-output").html()
			
			$("#tab-logs-error").css("display","none")
			$("#error_log").attr("class", "elm-header-log")
			
			$("#tab-logs-scipion").css("display","none")
			$("#scipion_log").attr("class", "elm-header-log")
			
		    break;
		    
		case "error_log":
			
			$("#tab-logs-output").css("display","none")
			$("#output_log").attr("class", "elm-header-log")
			
			$("#tab-logs-error").css("display","")
			$("#error_log").attr("class", "elm-header-log_selected")
			html = $("#tab-logs-error").html()
			
			$("#tab-logs-scipion").css("display","none")
			$("#scipion_log").attr("class", "elm-header-log")
		
		    break;
		    
		case "scipion_log":
			
			$("#tab-logs-output").css("display","none")
			$("#output_log").attr("class", "elm-header-log")
			
			$("#tab-logs-error").css("display","none")
			$("#error_log").attr("class", "elm-header-log")
			
			$("#tab-logs-scipion").css("display","")
			$("#scipion_log").attr("class", "elm-header-log_selected")
			html = $("#tab-logs-scipion").html()
		
		    break;
	}
	
	$("#externalTool").attr("href","javascript:customPopupFileHTML('" + html +"','"+ log_type + "',1024,768)")

}


function fillUL(type, list, ulId, icon) {
	/*
	 * Fill an UL element with items from a list items should contains id and name
	 * properties
	 */
	var ul = $("#" + ulId);
	ul.empty();
	for ( var i = 0; i < list.length; i++) {
		
//		inihtml = "<table><tr><td><strong>Attribute</strong></td><td><strong>&nbsp;&nbsp;&nbsp;Info</strong></td></tr>"
		
		// Visualize Object
		var visualize_html = '<a href="javascript:launchViewer(' + list[i].id
		+ ');"><i class="fa ' + icon + '" style="margin-right:10px;"></i>'+ list[i].name

		if(type=="input"){
			visualize_html += ' (from ' + list[i].nameId +')</a>'
		}
		else if(type=="output"){
			visualize_html += '</a>'
		}
		
		// Edit Object
		var edit_html = '<a href="javascript:editObject('+ list[i].id + ');"> '+
		'<i class="fa fa-pencil" style="margin-left:0px;"></i></a>'
		
//		endhtml = "</table>"
		
		ul.append('<li>' + visualize_html + edit_html +"&nbsp;&nbsp;&nbsp;" +list[i].info+ '</li>');
//		ul.append(inihtml+'<tr><td>' + visualize_html +"</td><td>&nbsp;&nbsp;&nbsp;" + list[i].info + edit_html + '</td></tr>' + endhtml);
		
	}
}

function updateButtons(id, mode){
	/*
	 * Function to update the buttons in the toolbar after choose a new protocol run.
	 */
	
	 if(mode == undefined){
		 mode="single";
	 }
	
	 switch (mode){
	 
	 	case "single":
	 		
			// Action Edit Button
			$("a#editTool").attr('href',
			'javascript:popup("/form/?protocolId=' + id + '")');
			$("span#editTool").show();
			
			// Action Copy Button
			$("a#copyTool").attr('href',
			'javascript:popup("/form/?protocolId=' + id + '&action=copy' + '")');
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
	 		
	 		var list_id = id
	 	
			// Action Edit Button
			$("a#editTool").attr('href','#');
			$("span#editTool").hide();
	 		
	 		// Action Copy Button
			$("a#copyTool").attr('href',
				'javascript:copyProtocol("' + list_id + '")');
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

function updateRow(id, elm, row){	
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
			dismarkSingleNodeGraph($(aux))
		}
		row.attr('value', id);
		$("div#graphActiv").attr("data-option", selected);
		markSingleNodeGraph(elm)
	}
}


function graphON(graph, icon_graph, list, icon_list){
	/*
	 * Function to disable the list view and enable the tree graph view.
	 */
	
	// Graph ON
	graph.attr("data-mode", "active");
//	graph.removeAttr("style");
	graph.show();
	graph.css("margin-top","-50em")
	icon_graph.hide();
	
	// Table OFF
	list.attr("data-mode", "inactive");
//	list.attr("style", "display:none;");
	list.hide();
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
//	list.removeAttr("style");
	list.show();
	icon_list.hide();
	
	// Graph OFF
	graph.attr("data-mode", "inactive");
//	graph.attr("style", "display:none;");
	graph.hide();
	$("img#loading").hide();
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
	
	changeStatusGraph(status, graph, graphTool, list, listTool);
		
	// Graph will be painted once
	if (graph.attr("data-time") == 'first') {
		callPaintGraph();
		graph.attr("data-time", "not");
	} 
	
	// Keep the consistency about the selected elements between
	// the list and graph views.
	transposeElmMarked(status);
	

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
		url : "/update_graph_view/?status=" + status,
		async: false
	});
}

function editObject(objId){
	$.ajax({
		type : "GET",
		url : '/get_attributes/?objId='+ objId,
		dataType: "text",
		success : function(text) {
			var res = text.split("_-_")
			label = res[0]
			comment = res[1]
			editObjParam(objId, "Label", label, "Comment", comment, "Describe your run here...","object")
		}
	});
}

function deleteProtocolForm(id, mode) {
	/*
	 * Dialog web form based in messi.js to verify the option to delete.
	 */
	
	 var msg = "</td><td class='content' value='"+ id
	 
	 switch(mode){
	 	case "single":
			msg += "'><strong>ALL DATA</strong> related to this <strong>protocol run</strong>"
			break;
	 	case "multiple":
			msg += "'><strong>ALL DATA</strong> related to the <strong>selected protocols run</strong>"
	 		break;
	 }
	 msg += " will be <strong>DELETED</strong>. Do you really want to continue?</td></tr></table>";
	 
	warningPopup('Confirm DELETE',msg, 'deleteProtocol')
	
}


function deleteProtocol(elm) {
	/*
	 * Method to execute a delete for a protocol
	 */
	var id = elm.attr('value');
	
	$.ajax({
		type : "GET",
		url : "/delete_protocol/?id=" + id,
		dataType : "json",
		async :false,
		success : function(json) {
			if(json.errors != undefined){
				// Show errors in the validation
				errorPopup('Errors found',json.errors);
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

function copyProtocol(id){
	$.ajax({
		type : "GET",
		url : "/copy_protocol/?id=" + id,
		dataType : "json",
		async :false,
		success : function(json) {
			if(json.errors != undefined){
				// Show errors in the validation
				errorPopup('Errors found',json.errors);
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
			
	warningPopup('Confirm STOP',msg, 'stopProtocol');
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

//** REFRESHING FUNCTIONALITIES *****************************************************/

function refreshRuns(mode){
	/*
	 * Method to update the run list/graph
	 */
	
	$(function() {
		$.ajax({
			url : '/run_table_graph/',
			success : function(data) {
			
				if (typeof data == 'string' || data instanceof String){
					
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
				else {
					for (var x=0;x< data.length;x++){
						var id = data[x][0];
						var status = data[x][2];
						var time = data[x][3];

						checkStatusNode(id, status)
						checkStatusRow(id, status, time)

					}
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

function checkStatusNode(id, status){
	node = $("#graph_"+ id);
	
	if(status.indexOf('running')==0){
		node.css("background-color", "#FCCE62");
	}
	if(status.indexOf('failed')==0){
		node.css("background-color", "#F5CCCB");
	}
	if(status.indexOf('finished')==0){
		node.css("background-color", "#D2F5CB");
	}
	node.find("#nodeStatus").html(status);
		
}

function checkStatusRow(id, status, time){
	row = $("tr#"+ id);
	
	row.find(".status").html(status);
	row.find(".time").html(time);
		
}

