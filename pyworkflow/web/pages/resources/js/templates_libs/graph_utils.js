/**
 * Graph methods and variables to use with jsPlumb plugin
 * 
 * callPaintGraph();
 * paintBox(nodeSource, id, msg);
 * addStatusBox(nodeSource, id, status);
 * connectNodes(elm1, elm2);
 * 
 **/

// Setting up drop options
var targetDropOptions = {
	tolerance : 'touch',
	hoverClass : 'dropHover',
	activeClass : 'dragActive'
};

// Setting up a Target endPoint
var targetColor = "black";
// var targetColor = "black";
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
// var sourceColor = "black";
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

function callPaintGraph() {
	// Draw the boxes
	var nodeSource = $("div#graphActiv");
	var status = "finished";
	var aux = [];

	// Paint the first node
	paintBox(nodeSource, "graph_PROJECT", "PROJECT", "");
	var width = $("div#" + "graph_PROJECT" + ".window").width();
	var height = $("div#" + "graph_PROJECT" + ".window").height();
	aux.push("PROJECT" + "-" + width + "-" + height);

	// Paint the other nodes (selected include)
	$("tr.runtr,tr.selected").each(function() {
		var id = jQuery(this).attr('id');
		var idNew = "graph_" + id;

		var name = jQuery(this).attr('data-name');

		paintBox(nodeSource, idNew, name, status);
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
				addStatusBox(nodeSource, "graph_" + json[i].id,
						json[i].status);
				$("div#graph_" + json[i].id + ".window").attr(
						"style",
						"top:" + top + "px;left:" + left
								+ "px;background-color:"
								+ json[i].color + ";");
			}
			// After all nodes are positioned, then create the edges
			// between
			// them

			for ( var i = 0; i < json.length; i++) {
				for ( var j = 0; j < json[i].childs.length; j++) {
					var source = $("div#graph_" + json[i].id
							+ ".window");
					var target = $("div#graph_" + json[i].childs[j]
							+ ".window");
					connectNodes(source, target);
				}
			}
			// If you choose first a element in the table, the
			// equivalent node
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

	if (id != "graph_PROJECT") {
		var objId = id.replace("graph_", "");
		var href = "javascript:popup('/form/?protocolId=" + objId + "')";
		var projName = $("div#graphActiv").attr("data-project");
		var onclick = "updateTabs('" + projName + "', '" + objId
				+ "',($(this)))";
		var aux = '<div class="window" style="" onclick="' + onclick + '" id="'
				+ id + '"><a href="' + href + '"><strong>' + msg
				+ '</strong></a><br /></div>';
	} else {
		var aux = '<div class="window" style="" id="' + id + '"><strong>' + msg
				+ '</strong><br />' + "" + '</div>';
	}

	nodeSource.append(aux);

	var oldSelect = $("div#graphActiv").attr("data-option");
}

function addStatusBox(nodeSource, id, status) {
	$("div#" + id + ".window").append(status);
}

function connectNodes(elm1, elm2) {
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

// function putEndPoints(elm) {
// jsPlumb.addEndpoint($(elm + ".window"), {
// anchor : "TopCenter"
// }, targetEndpoint);
// jsPlumb.addEndpoint($(elm + ".window"), {
// anchor : "BottomCenter"
// }, sourceEndpoint);
// }