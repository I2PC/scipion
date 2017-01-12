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

 /** METHODS ******************************************************************/
	
function launchToolbarObjTree(id, elm, multi) {
	/*
	 * Toolbar used in the data content template for graph view
	 */
	
	if (multi){
		if (elm.attr("selected") == "selected"){
			dismarkSingleNodeGraph(elm);
			$("div#graphActiv").removeAttr("data-option");
			updateObjButtons(id, "multiple")
		} else {
			markSingleNodeGraph(elm);
			$("div#graphActiv").attr("data-option", elm.attr("id"));
		}
		updateButtons(id, "multiple")
	} else {
		disableMultipleMarkGraph(id);
		updateTree(id, elm, $("div#toolbar"));
		// Update the content for the tabs
		updateObjButtons(id, "single")
		updateObjTabs(id);
		
	}
}

function updateObjTabs(id) {
	/*
	 * Fill the content of the summary tab for a object selected 
	 */
	var URL = getSubDomainURL() + '/object_info/?objectId=' + id
	
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		success : function(json) {

			// Edit Object
			var edit_html = '<a href="javascript:editObject('+ id + ');"> '+
			'<i class="fa fa-pencil" style="margin-left:0px;"></i></a>'
			
			// Label
			var label = $("#obj_label");
			label.empty();
			var value_label = $("#graph_"+id).attr("data-label")
			label.append(value_label + "&nbsp;&nbsp;&nbsp;" + edit_html);
			
			// Info
			var info = $("#obj_info");
			info.empty();
			info.append(json.info);
			
			// Created
			var created = $("#obj_created");
			created.empty();
			created.append(json.created);
			
			//Comment
			var comment = $("#obj_comment");
			comment.empty();
			comment.append(json.comment);
		}
	});
}

function updateObjButtons(id, mode){
	/*
	 * Function to update the buttons in the toolbar after choose a new protocol run.
	 */	
	 switch (mode){
	 	case "single":
			// Action Edit Button
			$("a#editTool").attr('href', 'javascript:editObject('+ id + ');');
			$("span#editTool").show();
			
			// Show toolbar
			$("div#toolbar").show(); 
			break;
			
	 	case "multiple":
	 		// Hide toolbar
	 		$("div#toolbar").hide(); 
	 		break;
	 }
}


function callPaintObjTree(){
	var URL = getSubDomainURL() + '/object_tree/'
	$.ajax({
		type : "GET",
		url : URL,
		dataType: "html",
		async: false,
		success : function(data) {
			$("ul#browser").html(data)
		}
	});
	$("#browser").treeview()
}

