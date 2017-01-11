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
 * Graph methods and variables to use with jsPlumb plugin
 * 
 * ATTRIBUTES LIST:
 * 
 * var targetDropOptions
 * var targetColor
 * var targetEndpoint
 * var sourceColor
 * var sourceEndpoint
 * 
 * METHODS LIST:
 * 
 * function callPaintGraph()
 * 	->	This function paint the protocol graph in the template project_content.
 * 
 * function paintBox(node_source, id_node, text_node)
 * 	->	Function to paint a box like a node inside the protocol graph.
 * 	    The node source is passed by arguments.
 *			usage example: paintBox(nodeSource, "protocol_new", "New Protocol");
 * 
 * function addStatusBox(id_node, status_text)
 * 	->	Function add a new status in a node from the protocol graph.
 * 			usage example: addStatusBox("protocol_new", "finished"); 
 * 
 * function connectNodes(node_1, node_2)
 * 	->	Function to connect two nodes with a line using jsPlumb.
 * 
 ******************************************************************************/
 
/** ATTRIBUTES ****************************************************************/
	
// Setting up drop options
var targetDropOptions = {
	tolerance : 'touch',
	hoverClass : 'dropHover',
	activeClass : 'dragActive'
};

// Setting up a Target endPoint
var targetColor = "black";

var targetEndpoint = {
	endpoint : [ "Dot", {
		radius : 5
	} ],
	paintStyle : {
		fillStyle : targetColor
	},
	// isSource:true,
	scope : "green dot",
	connectorStyle : {
		strokeStyle : targetColor,
		lineWidth : 2
	},
	connector : [ "Bezier", {
		curviness : 2
	} ],
	maxConnections : 100,
	isTarget : true,
	dropOptions : targetDropOptions
};

// Setting up a Source endPoint
var sourceColor = "black";

var sourceEndpoint = {
	endpoint : [ "Dot", {
		radius : 5
	} ],
	paintStyle : {
		fillStyle : sourceColor
	},
	isSource : true,
	scope : "green dot",
	connectorStyle : {
		strokeStyle : sourceColor,
		lineWidth : 2
	},
	connector : [ "Bezier", {
		curviness : 2
	} ],
	maxConnections : 100
	// isTarget:true,
	// dropOptions : targetDropOptions
};

/** METHODS *******************************************************************/

function callPaintGraph(nodeSource, type) {
	/*
	 * This function paint the protocol graph in the template project_content.html
	 */ 
	
	var nameId = "graph_"
	if (type == "small"){
		nameId = "graphSmall_"
	}
	
	var status = "finished";
	var aux = [];

	// Paint the first node
	paintBox(nodeSource, nameId, nameId + "PROJECT", "PROJECT", "runs");
	var width = $("div#" + nameId + "PROJECT" + ".window").width();
	var height = $("div#" + nameId + "PROJECT" + ".window").height();
	aux.push("PROJECT" + "-" + width + "-" + height);

	// Paint the other nodes (selected include)
	$("tr.runtr,tr.selected").each(function() {
		var id = jQuery(this).attr('id');
		var idNew = nameId + id;

		var name = jQuery(this).attr('data-label');
		if (name==""){
			name = jQuery(this).attr('data-name');
		}

		if(type == "normal"){
			paintBox(nodeSource, nameId, idNew, name, "runs");
		} else if (type == "small"){
			paintBox(nodeSource, nameId, idNew, id, "runs");
		}
		
		var width = $("div#" + idNew + ".window").width();
		var height = $("div#" + idNew + ".window").height();

		aux.push(id + "-" + width + "-" + height);
	});	
	
	var URL = getSubDomainURL() + '/project_graph/?list=' + aux
	// Move and connect the graph nodes
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		async: false,
		success : function(json) {
			processNodes(json, "runs", nameId, type)
		}
	});
	jsPlumb.draggable($(".window"));
}


function callPaintObjGraph(){
	/*
	 * This function paint the object graph in the template data_content.html
	 */ 
	var nodeSource = $("div#graphActiv");
	var URL = getSubDomainURL() + '/elements_graph/'
	// Get the objects information for be painted
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		async: false,
		success : function(json) {
			// Paint the nodes
			var aux = [];
			$.each(json, function(i, item) {
				var id = "graph_" + item.id
				var label = item.label
				
				// Draw box
				paintBox(nodeSource, "graph_", id, label, "objects");
				
				if(item.id!= 'PROJECT'){
					// Get the size
					var width = $("div#" + id + ".window").width();
					var height = $("div#" + id + ".window").height();
					aux.push(item.id + "-" + width + "-" + height);
				}
			});
			
			// Move and connect the graph nodes
			var URL = getSubDomainURL() + '/object_graph/?list=' + aux
			$.ajax({
				type : "GET",
				url : URL,
				dataType : "json",
				async: false,
				success : function(json) {
					processNodes(json, "objects", "graph_")
				}
			});
		}
	});
	jsPlumb.draggable($(".window"));
}


function paintBox(nodeSource, nameId, id, msg, mode) {
	/*
	 * Function to paint a box like a node inside the protocol graph.
	 * The node source is passed by arguments.
	 * 
	 * FIX 1: Added the css property (display:none) to the divs to not show
	 * in first instance the boxes in bad position, before to be processed.
	 * 
	 */
	var auxId = nameId + "PROJECT";
	if (id != auxId) {
		var objId = id.replace(nameId, "");
		
		switch(mode){
			case "runs":
				var href = "javascript:customPopup('/form/?protocolId=" + objId + "',620,591)";
				var onclick = "launchToolbarTree('" + objId	+ "', $(this), isCtrlPress(event))";
				break;
			case "objects":
				var href = "javascript:launchViewer(" + objId +")";
				var onclick = "launchToolbarObjTree('" + objId	+ "', $(this), isCtrlPress(event))";
		}
		
		var projName = nodeSource.attr("data-project");
		var aux = '<div class="window" style="display:none;" onclick="' + onclick + '" id="'
				+ id + '" data-label="'+ msg +'"><a href="' + href + '"><strong>' + msg
				+ '</strong></a><br/><span id="nodeStatus" data-val=""></span></div>';	
	} else {
		var aux = '<div class="window" style="display:none;" id="' + id + '"><strong>' + msg
				+ '</strong><br /></div>';
	}
	nodeSource.append(aux);
}


function processNodes(json, mode, nameId, type){
	// Iterate over the nodes and position in the screen
	// coordinates should come in the json response
	
	positionNodes(json, mode, nameId, type)
	
	// After all nodes are positioned, then create the edges
	// between them
	
	putEdges(json, nameId)
}

function positionNodes(json, mode, nameId, type){
	for ( var i = 0; i < json.length; i++) {
		var top = json[i].y*0.8;
		var left = json[i].x;
		
		switch(mode){
			case "objects":
				var style = "top:" + top * 1.1
						+ "px;left:" + left*1.2
						+ "px;background-color:"+ json[i].color + ";"
						+ "padding:8px;"
				break;
				
			case "runs":
				if (type!="small"){
					addStatusBox(nameId + json[i].id,	json[i].status);
				}
						
				var style = "top:" + top
					+ "px;left:" + left
					+ "px;background-color:"+ json[i].color + ";"
				break;
		}
		
		$("div#"+ nameId + json[i].id + ".window").attr(
				"style", style);
	}
}

function addStatusBox(id, status) {
	/*
	 * Function add a new status in a node from the protocol graph.
	 */
	$("div#" + id + ".window").find("#nodeStatus").html(status);
}

function putEdges(json, nameId){
	for ( var i = 0; i < json.length; i++) {
		for ( var j = 0; j < json[i].childs.length; j++) {
			var source = $("div#"+ nameId + json[i].id
					+ ".window");
			var target = $("div#"+ nameId + json[i].childs[j]
					+ ".window");
			connectNodes(source, target);
		}
	}
}

function connectNodes(elm1, elm2) {
	/*
	 * Function to connect two nodes with a line using jsPlumb.
	 */
	
	// alert($(elm1).attr('id') + " - " + $(elm2).attr('id'));
	var a = jsPlumb.addEndpoint($(elm1), {
		anchor : "Center"
	}, sourceEndpoint);
	var b = jsPlumb.addEndpoint($(elm2), {
		anchor : "Center"
	}, targetEndpoint);

	jsPlumb.connect({
		source : a,
		target : b
	});
}
