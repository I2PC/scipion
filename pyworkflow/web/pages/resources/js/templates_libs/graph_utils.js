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

function callPaintGraph() {
	/*
	 * This function paint the protocol graph in the template project_content.html
	 */ 
	
	// Draw the boxes
	var nodeSource = $("div#graphActiv");
	var status = "finished";
	var aux = [];

	// Paint the first node
	paintBox(nodeSource, "graph_PROJECT", "PROJECT");
	var width = $("div#" + "graph_PROJECT" + ".window").width();
	var height = $("div#" + "graph_PROJECT" + ".window").height();
	aux.push("PROJECT" + "-" + width + "-" + height);

	// Paint the other nodes (selected include)
	$("tr.runtr,tr.selected").each(function() {
		var id = jQuery(this).attr('id');
		var idNew = "graph_" + id;

		var name = jQuery(this).attr('data-name');

		paintBox(nodeSource, idNew, name);
		var width = $("div#" + idNew + ".window").width();
		var height = $("div#" + idNew + ".window").height();

		aux.push(id + "-" + width + "-" + height);
	});

	$.ajax({
		type : "GET",
		url : '/project_graph/?list=' + aux,
		dataType : "json",
		success : function(json) {
		
			// Iterate over the nodes and position in the screen
			// coordinates should come in the json response
			
			for ( var i = 0; i < json.length; i++) {
				var top = json[i].y*0.8;
				var left = json[i].x;
				addStatusBox("graph_" + json[i].id,	json[i].status);
				$("div#graph_" + json[i].id + ".window").attr(
						"style",
						"top:" + top + "px;left:" + left
								+ "px;background-color:"
								+ json[i].color + ";");
			}
			
			// After all nodes are positioned, then create the edges
			// between them

			for ( var i = 0; i < json.length; i++) {
				for ( var j = 0; j < json[i].childs.length; j++) {
					var source = $("div#graph_" + json[i].id
							+ ".window");
					var target = $("div#graph_" + json[i].childs[j]
							+ ".window");
					connectNodes(source, target);
				}
			}
			
			// If you choose first a element in the table, the equivalent node
			// must be flashlighted in the graph
			
			if ($("tr.selected").attr("id") != undefined) {
				var selected = "graph_" + $("tr.selected").attr("id");
				$("div#graphActiv").attr("data-option", selected);
				var elm = $("div#" + selected + ".window");
				var aux = elm.attr("style");
				aux += "border:2.5px solid Firebrick;"
				elm.attr("style", aux);
			}
		}
	});
	jsPlumb.draggable($(".window"));
}

function paintBox(nodeSource, id, msg) {
	/*
	 * Function to paint a box like a node inside the protocol graph.
	 * The node source is passed by arguments.
	 */

	if (id != "graph_PROJECT") {
		var objId = id.replace("graph_", "");
		var href = "javascript:popup('/form/?protocolId=" + objId + "')";
		var projName = $("div#graphActiv").attr("data-project");
//		var onclick = "updateTabs('" + projName + "', '" + objId + "',($(this)))";
		var onclick = "launchToolbarTree('" + objId	+ "',($(this)))";
		var aux = '<div class="window" style="" onclick="' + onclick + '" id="'
				+ id + '"><a href="' + href + '"><strong>' + msg
				+ '</strong></a><br/><span id="nodeStatus" data-val="hola"></span></div>';	
	} else {
		var aux = '<div class="window" style="" id="' + id + '"><strong>' + msg
				+ '</strong><br /></div>';
	}
	// + '</strong><br /><span id="nodeStatus" data-val="hola"></span></div>';
	nodeSource.append(aux);
//	var oldSelect = $("div#graphActiv").attr("data-option");
}

function addStatusBox(id, status) {
	/*
	 * Function add a new status in a node from the protocol graph.
	 */
	
//	$("div#" + id + ".window").append(status);
	$("div#" + id + ".window").find("#nodeStatus").html(status);
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
