/*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
 * 			   
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
 * Methods used in the 'showj_menu' and 'showj_bottom_menu' templates. 
 * 
 * METHODS LIST:
 * 
 * function getListEnabledItems(mode)
 * 	-> 	Obtains the list of the items selected like enable.
 * 
 * function getListSelectedItems(mode)
 * 	-> 	Obtains the list of the items selected in the template.
 * 
 * function changeMode(modeNew)
 * 	->	Function to change the visualization mode.
 * 
 * function markSelectedItems(mode, list)
 * 	->	Function to mark the items what the list contains.
 * 
 * function updateEnabledItems(mode, list)
 * 	->	Function to update the enabled items what the list contains.
 * 
 * function updateCheckboxDataTable(element, mode)
 * 	->	This method was done to keep updated the data core in
 *		the datatable library used to show in the table mode.
 *
 * function updateListSession(id, attr)
 * 	->	Method to update the session variable about list depending on attr.
 * 
 * function updateList(id, list)
 * 	->	Method to add/remove the element identified into the list.
 * 
 * function updateSelectedTemplate(list)
 * 	->	Method to update the template with the selected elements.
 * 
 * function updateEnabledTemplate(list, persist)
 * 	->	Method to update the template with the enables elements.
 * 		Option to persist the element into the oTable object available.
 * 
 ******************************************************************************/
	 
function getListEnabledItems(mode){
	/*
	 * Obtains the list of the items selected like enable.
	 */
	var list_enabled=[];
    switch(mode){
        case "gallery":
            // To get the enabled images
            $.each($("img.enabledGallery"), function(){
                var display = $(this).css("display");
                if (display=='block'){
                	var id = $(this).attr('id').split("___").pop();
                	list_enabled.push(id);
                }
            });
            break;
            
        case "table":
            // To get the enabled rows
            $.each($("input[id^='enabled___']:not(:checked)"), function(){
           		var id = $(this).attr("id").split("___").pop();
           		list_enabled.push(id);
            });
            break;
    }
    return list_enabled
}

function getListSelectedItems(mode){
	/*
	 * Obtains the list of the items selected in the template.
	 */
	var list_selected=[];
    switch(mode){
        case "gallery":
            // To get the selected images
            $.each($("div.image_selected"), function(){
                var id = $(this).attr("id").split("___").pop()
                list_selected.push(id)
            });
            break;
            
        case "table":
            // To get the selected rows
            $.each($("tr.row_selected"), function(){
                var id = $(this).attr("id").split("___").pop()
                list_selected.push(id)
            });
            break;
    }
    return list_selected
}

function changeMode(modeNew){
	/*
	 * Function to change the visualization mode.
	 */
	var form = document.forms['showjForm']
	var modeOld = form.mode.value
	 
	if (modeNew != modeOld){
		
		// Selected Items
		var listSelected = $("input#listSelected").val();
		if(listSelected == 0 || listSelected == '0'){
			// First execution
			listSelected = getListSelectedItems(modeOld);
		}
		
		// Enabled Items
		var listEnabled = $("input#listEnabledItems").val();
		if(listEnabled == 0 || listEnabled == '0'){
			// First execution
			listEnabled = getListEnabledItems(modeOld);
		}
		form.listEnabledItems.value = listEnabled;
		
		
		// Submit form
		form.mode.value = modeNew;
 		form.submit();
	}
}

function markSelectedItems(mode, list){
	/*
	 * Function to mark the items what the list contains.
	 */
	switch(mode){
	    case "gallery":
	        // Came from the table mode
	        for (var x=0;list.length>x;x++){
	        	$("div#img_container___"+ list[x]).addClass("image_selected");
	        }
	        break;
	        
	    case "table":
	        // Came from the gallery mode
	    	updateSelectedTemplate(list)
	        break;
	}
}

function updateEnabledItems(mode, list){
	/*
	 * Function to update the enabled items what the list contains.
	 */
	switch(mode){
	    case "gallery":
	        // Came from the table mode
	        for (var x=0;list.length>x;x++){
	        	var elm = $("img#enabled___"+ list[x]);
	        	elm.css("display","inline");
	        	elm.addClass("selected");
	        }
	        break;
	        
	    case "table":
	        // Came from the gallery mode
	    	updateEnabledTemplate(list, true)
	        break;
	}
}

function updateCheckboxDataTable(element, mode){
	/*
	 * This method was done to keep updated the data core in
	 * the datatable library used to show in the table mode. 
	 */
	var aPos = oTable.fnGetPosition(element.parents('td').get(0))
	
//	  @param {object|array|string} mData Data to update the cell/row with
//    @param {node|int} mRow TR element you want to update or the aoData index
//    @param {int} [iColumn] The column to update (not used of mData is an array or object)
//    @param {bool} [bRedraw=true] Redraw the table or not
//    @param {bool} [bAction=true] Perform pre-draw actions or not
//    @returns {int} 0 on success, 1 on error
	oTable.fnUpdate(mode, aPos[0], 1, true);
	
}

function updateListSession(id, attr){
	/*
	 * Method to update the session variable about list depending on attr.
	 */
	switch(attr){
		case "enabled":
			updateList(id, "listEnabledItems");
			break;
		
		case "selected":
			updateList(id, "listSelectedItems");
			break;
	}
}


function updateList(id, list){
	/*
	 * Method to add/remove the element identified into the list.
	 */
	var id = id.split("___").pop();
	var listItems = $("#" + list).val().split(",")
	
	var enc=-1;
	
	var x=0;
	do{
		if (listItems[x] == id){
			enc=x;
		}
		x++;
	} while(x<listItems.length && enc==-1);
	
	
	if(enc!=-1){
		listItems.splice(enc, 1);
	}else{
		listItems.push(id)
	}
	
	if(listItems[0]==""){
		listItems.shift()
	}
	
	$("#" + list).val(listItems)
}

function updateSelectedTemplate(list){
	/*
	 * Method to update the template with the selected elements. 
	 */
	for (var x=0;list.length>x;x++){
		$("tr#row_container___"+ list[x]).addClass("row_selected");
	}
}

function updateEnabledTemplate(list, persist){
	/*
	 * Method to update the template with the enables elements.
	 * Option to persist the element into the oTable object available.
	 */
	if(persist == undefined)
		persist = false;
	 
	for (var x=0;list.length>x;x++){
		var elem = $("input#enabled___"+ list[x])
		elem.prop("checked", false);
		
		if(persist){
			//Fix to keep the datatable updated
    		updateCheckboxDataTable(elem, "False");
		}
		
	}
	
}


function createSubset(){
	var listSubset = $("#listEnabledItems").val()
	
	alert("To do a subset without the elements:" + listSubset);
	
}

