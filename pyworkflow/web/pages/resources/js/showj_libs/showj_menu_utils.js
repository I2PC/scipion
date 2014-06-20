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
		var listSeletected = getListSelectedItems(modeOld);
		var listEnabled = getListEnabledItems(modeOld);
		form.listSelectedItems.value = getListSelectedItems(modeOld);
		form.listEnabledItems.value = getListEnabledItems(modeOld);
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
	    	for (var x=0;list.length>x;x++){
	    		$("tr#row_container___"+ list[x]).addClass("row_selected");
	    	}
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
	        	$("img#enabled___"+ list[x]).css("display","block");
	        }
	        break;
	        
	    case "table":
	        // Came from the gallery mode
	    	for (var x=0;list.length>x;x++){
	    		$("input#enabled___"+ list[x]).prop("checked", false);
	    	}
	        break;
	}
}

function createSubset(){
	var mode = document.forms['showjForm'].mode.value
	var listSubset = getListEnabledItems(mode)
	
	alert("To do a subset with the elements:" + listSubset);
	
}