 /*****************************************************************************
 *
 * Authors:    Adrian Quintana (aquintana@cnb.csic.es)
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
 * Methods used in the showj_table template 
 * 
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * 
 * function renderElements(nRow, aData)
 * -> Render elements on table (Image and checkbox). This method is called from datatable script
 * 	  nRow is the row to be displayed
 *    aData the row content
 *     
 * function initializeColumnHeader()
 * -> Initialize icons in column header (render image, disable column and set column editable)
 * 
 * function getHeaderWithIcon(text, columnLayoutProperties) 
 * -> Generate html for renderable, editable and disable icon for each column header
 *    text is the text to be displayed in the header
 *    columnLayoutProperties defines the icons for column header
 *    
 * function getColumnsDefinition()  
 * -> Set the columns property for datatable script (visible, title, target...)
 * 
 * function enableDisableColumn(event, element)
 * -> Enable/disable column 
 *    Disable datatable column when event triggered from icon disable column header
 *    event is used to stop propagation
 *    element is used to get the column number 
 * ATTENTION:Currently this is not being used as delete icon has been removed
 * 
 * function enableDisableRenderColumn(event, element)
 * -> Enable/disable render property to column
 *    When column cell is image, it toggles between img element and path text
 *    event is used to stop propagation
 *    element is used to get the column number
 *    
 * function enableDisableEditableColumn(event, element)  
 * -> Enable/disable editable property to column
 *    When column cell is text or a number, it toggles between text element and input textfield
 *    event is used to stop propagation
 *    element is used to get the column number
 * 
 * function setElementsEditable(elements)
 * -> Make all the elements editable. Onclick elements will change into textfield (jeditable)
 *    If no elements are given, all dataset is configured
 *    
 * function valueChange(element)
 * -> Keep changed element in a global variable   
 * 
 * function initializeTableWidth()
 * -> If table content is wider than web width, set table content to page width.
 * 	  Executed on load in order to make scroll bars visible
 * 
 * function initializeGoToEvent()
 * -> Initialize events for Goto textfield
 * 	  Goto textfield allows to select a row in the table
 * 
 * function updateRowSelection()
 * -> Select the row set in the Goto textfield and scroll to the row
 * 
 * function initializeSelectionRowEvent()
 * -> Initialize onclick events for row.
 *    Single selection and special key controls (ctrl, shift, single click)
 * 
 * function initializeMultipleSelectionTool()
 * -> Initialize hover events for row selection.
 *    It will display the Multiple selection tool
 *    
 * function showTableConfig()
 * -> Display and initialize div for configuring table layout
 * 
 * function saveTableConfiguration()
 * -> Save new table configuration and redraw data table
 * 
 * function multipleEnableDisableImage(mode)
 * -> Enables/disable a set of rows.
 *    Mode can be enable or disable
 *    
 * function multipleSelect(mode)
 * -> Select a set of images
 *    Mode can be all, from (i.e. from the image that triggered the event until the end)
 *    or to (i.e. from the beginning until the image that triggered the event)
 *       
 ******************************************************************************/

 /** METHODS ******************************************************************/

function renderElements(nRow, aData){
	 columnId = 0
	 invisibleColumns = 0 
	 for (var i in jsonTableLayoutConfiguration.columnsLayout){ 
		 var columnLayoutConfiguration = jsonTableLayoutConfiguration.columnsLayout[i]
		 
		 if (columnLayoutConfiguration.columnLayoutProperties.visible){
	   		 if (typeof oTable!= 'undefined'){
	   			 columnId=oTable.fnGetColumnIndex(i)
	   			 columnIdReal=oTable.fnColumnIndexToVisible(columnId)
	   		 }
	   		 else{
	   			 columnIdReal=columnId-invisibleColumns
	   		 }
	   		 
	   		if (columnIdReal!=null){ 
	   			if (columnLayoutConfiguration.typeOfColumn == 'checkbox'){
		   			 checkbox_element = '<input type=\"checkbox\" onclick=\"valueChange(this);\" id=\"'+i+'___'+aData[0]+'\" '
		    			 if (aData[columnId] == 1){checkbox_element += 'checked>'}
		    			 else{checkbox_element += '>'}
		   			 $('td:eq('+columnIdReal+')', nRow).html(checkbox_element);
		   		}
	   			else if (columnLayoutConfiguration.columnLayoutProperties.renderable){
		   			$('td:eq('+columnIdReal+')', nRow).html(
		   					'<span style="display:none">'+aData[columnId]+'</span>'+
		   					'<img class=\"tableImages\" id=\"'+i+'___'+aData[0]+'\" src=\"/render_column/?renderFunc='+columnLayoutConfiguration.columnLayoutProperties.renderFunc+'&'+columnLayoutConfiguration.columnLayoutProperties.extraRenderFunc+'&image='+aData[columnId]+'\"/>' );
		   		 
	   			}
	   			else if (columnLayoutConfiguration.typeOfColumn == 'image'){
			   			$('td:eq('+columnIdReal+')', nRow).html( '<span>'+aData[columnId]+'</span>'+
			   					'<img style="display:none" class=\"tableImages\" id=\"'+i+'___'+aData[0]+'\" data-real_src=\"/render_column/?renderFunc='+columnLayoutConfiguration.columnLayoutProperties.renderFunc+'&'+columnLayoutConfiguration.columnLayoutProperties.extraRenderFunc+'&image='+aData[columnId]+'\"/>' );
			   			/* $('td:eq('+columnIdReal+')', nRow).html( '<span>'+aData[columnId]+'</span>' ); */
		   		 }
	   		}	
		 }
		 else{
			 invisibleColumns++;
		 }
		 columnId++;
	 }
}

function initializeColumnHeader(){
	headerRow=$("#data_table thead tr")
	var cols = oTable.fnSettings().aoColumns;
    for ( var x=0, xLen=cols.length ; x<xLen ; x++ )
    {
        var cellIndex = oTable.fnColumnIndexToVisible(x)
        if (cellIndex !=null){
        	$('th:eq('+cellIndex+')',headerRow).find("div").append(getHeaderWithIcon(cols[x].sTitle,jsonTableLayoutConfiguration.columnsLayout[cols[x].sTitle].columnLayoutProperties))
        	$('th:eq('+cellIndex+')',headerRow).attr("id",cols[x].sTitle+"___column_header")
        } 
    }
	
} 

function getHeaderWithIcon(text, columnLayoutProperties){
	var iconElements = '' 
//	 NAPA DE LUXE HABRIA QUE PONERLE UN MARGIN AL LABEL DEL HEADER
//	var iconElements = '<span style=\"margin:8px\">'+text+'</span>'
// 	if (columnLayoutProperties.visible && columnLayoutProperties.allowSetVisible){
//		iconElements+="<span class=\"css_right\"><img src=\"/resources/showj/delete.png\" class=\"arrowImage\" onclick=\"enableDisableColumn(event,this)\" ></span>"
//	} 
	
	if (columnLayoutProperties.allowSetEditable){
		iconElements+="<span class=\"css_right fa-stack\"><a id=\""+ text+"_editable_icon\" href='#' onclick=\"enableDisableEditableColumn(event,this);\"><i class=\"fa fa-pencil fa-stack-1x"
		iconElements+="\"></i><i id=\"banElement\" class=\"fa fa-stack-1x "
		if (columnLayoutProperties.renderable){iconElements+="fa-times"}
		iconElements+="\"></i></a></span>"
	}
	
	if (columnLayoutProperties.allowSetRenderable){
		iconElements+="<span class=\"css_right\"><a id=\""+ text+"_renderable_icon\" href='#' onclick=\"javascript:enableDisableRenderColumn(event,this);\"><i class=\"fa "
		if (columnLayoutProperties.renderable){iconElements+="fa-eye"}else{iconElements+="fa-eye-slash"}
		iconElements+="\"></i></a></span>"
	}
	return iconElements; 
}

function getColumnsDefinition(){
	jsonColumnsLayout=jsonTableLayoutConfiguration.columnsLayout
	columnId=0;
	var dataForTable = []
	for (var i in jsonColumnsLayout){ 
		var dataRowForTable = []
 		var columnLayoutConfiguration = jsonColumnsLayout[i]
 		if(!columnLayoutConfiguration.columnLayoutProperties.visible){
			dataRowForTable["bVisible"] = false
		} 
		dataRowForTable["sTitle"]= i
		dataRowForTable["aTargets"]=[columnId]
		dataForTable.push(dataRowForTable)
		columnId++;
	}
	
	return dataForTable;
}
 
function enableDisableColumn(event, element){
	//From the image element we get the column header index
	var thCell= $(element).closest('th')
	var iCol = $('th').index(thCell)
	var iCol2 = oTable.fnVisibleToColumnIndex(iCol)
	
	var bVis = oTable.fnSettings().aoColumns[iCol2].bVisible;
	oTable.fnSetColumnVis( iCol2, bVis ? false : true );

//	Update table layout configuration model
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].columnLayoutProperties.visible=!bVis
	
	//This will avoid column sorting
	event.stopPropagation() 
}
 
function enableDisableRenderColumn(event, element){
//	Switch button from on to off or viceversa
	$(element).find("i").toggleClass("fa-eye").toggleClass("fa-eye-slash")
	
	//From the image element we get the column header index
	var thCell= $(element).closest('th')
	var iCol = $('th').index(thCell)
	
	//We enable/disable render event
	var nTrs = oTable.fnGetNodes();
	$('td:nth-child('+(iCol+1)+')', nTrs).each(function(){
		if ($(element).find("i").hasClass("fa-eye") && $(this).find("img").attr('src')== undefined){
			/* $(this).find("img").attr("src",$(this).find("img").data('real_src')) */
			setImageSize(false)
			initializeImageLoad(false)
		}
		$(this).find("img").toggle()
		$(this).find("span").toggle()
	})

//	Update table layout configuration model
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].columnLayoutProperties.renderable=$(element).find("i").hasClass("fa-eye")
	
	//This will avoid column sorting
	event.stopPropagation() 
}
 
function enableDisableEditableColumn(event, element){
	$(element).find("#banElement").toggleClass("fa-times")
	
	//From the image element we get the column header index
	var thCell= $(element).closest('th')
	var iCol = $('th').index(thCell)
	
	//We enable/disable render event
	var nTrs = oTable.fnGetNodes();
	if (!$(element).find("#banElement").hasClass("fa-times")){
		$('td:nth-child('+(iCol+1)+')', nTrs).addClass('editable')
		setElementsEditable('td:nth-child('+(iCol+1)+')')
	}
	else{
		$('td:nth-child('+(iCol+1)+')', nTrs).removeClass('editable')
		$('td:nth-child('+(iCol+1)+')', nTrs).unbind('click');
	}
	$('td:nth-child('+(iCol+1)+')', nTrs).each(function(){
		$(this).toggleClass("editable")
	})
	
//	Update table layout configuration model
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].columnLayoutProperties.editable=!$(element).find("#banElement").hasClass("fa-times")
	
	//This will avoid column sorting
	event.stopPropagation()
	
}

function setElementsEditable(elements){
	if (elements==null){
		var nTrs = oTable.fnGetNodes();
		for (var label in jsonTableLayoutConfiguration.columnsLayout){ 
			columnId=oTable.fnGetColumnIndex(label)
   			columnIdReal=oTable.fnColumnIndexToVisible(columnId)
   			if (columnIdReal!=null){
				if (jsonTableLayoutConfiguration.columnsLayout[label].columnLayoutProperties.editable){
	   				setElementsEditable('td:nth-child('+(columnIdReal+1)+')')
				}
				else{
	   				$('td:nth-child('+(columnIdReal+1)+')', nTrs).unbind('click');
	   			}
   			}	
		}
	}
	else{
		$(elements, oTable.fnGetNodes()).editable(
			function(value, settings){ 
//				Get position and real column
				var aPos = oTable.fnGetPosition( this )
				columnIdReal=oTable.fnVisibleToColumnIndex(aPos[1])
				
//				Get label from column
				var label = oTable.fnSettings().aoColumns[columnIdReal].sTitle
				
//				NAPA DE LUX TENEMOS QUE PONERLE UN ID A LAS FILAS
				var aoData = oTable.fnSettings().aoData;
				var nTd = aoData[ aPos[0] ]._anHidden[ oTable.fnGetColumnIndex("id") ];
				var rowId =$(nTd).text()
				
//				Keep changes in global variable 
				if (!$("#saveButton").hasClass("buttonGreyHovered")){
					$("#saveButton").toggleClass("buttonGreyHovered")
				}
					
				changes[label+"___"+rowId]=value
			    return(value);
			  },
			{
			 "callback": function( sValue, y ) {
//				Update datatable model
				var aPos = oTable.fnGetPosition( this );
				columnIdReal=oTable.fnVisibleToColumnIndex(aPos[1])
				oTable.fnUpdate( sValue, aPos[0], columnIdReal );
			},
			"height": "14px"
		} );
	}	
}

function valueChange(element){
	element_value = ""
	if ($(element).is("input:checkbox")){element_value = $(element).is(":checked")}
	else{element_value = $(element).val()}
//	Keep changes in global variable
	changes[$(element).attr("id")]=element_value
}

function initializeTableWidth(){
	var tableWidth = $("section").width()
	if ($("#table_container").width()>tableWidth){
		$("#table_container").width(tableWidth);
	}	
}

function initializeGoToEvent(){
	$('#id_goto').click(function(){
		updateRowSelection();	
	});
	
	$('#id_goto').change(function(){
		updateRowSelection();
	});
}
function updateRowSelection(){
	oTable.fnDisplayRowWithIndex($("#id_goto").val()-1)
	var nodes = oTable.fnGetNodes();
	selectedElement=nodes[$("#id_goto").val()-1]
	selectedElement.click();
	
	$('html, body').animate({
        scrollTop: $(selectedElement).offset().top
    }, 2000);
}

var lastChecked;
function initializeSelectionRowEvent(){
	$('#data_table').on("click","tr", function(event) {

		var start = $('#data_table tbody tr').index(this);
		if (start>=0){
		    if(!lastChecked) {
		        lastChecked = this;
		    }
		     
		    if(event.shiftKey) {
		        var end = $('#data_table tbody tr').index(lastChecked);
		 
		        for(i=Math.min(start,end);i<=Math.max(start,end);i++) {
		            if (!$('#data_table tbody tr').eq(i).hasClass('row_selected')){
		                $('#data_table tbody tr').eq(i).addClass("row_selected");
		            }
		        }
		         
		        // Clear browser text selection mask
		        if (window.getSelection) {
		            if (window.getSelection().empty) {  // Chrome
		                window.getSelection().empty();
		            } else if (window.getSelection().removeAllRanges) {  // Firefox
		                window.getSelection().removeAllRanges();
		            }
		        } else if (document.selection) {  // IE?
		            document.selection.empty();
		        }
		    } else {
//		    	$(lastChecked).removeClass('row_selected');
		    	if (!event.metaKey && !event.ctrlKey){
		    		$('tr').each(function (){
						$(this).removeClass('row_selected');
					});
			    }
	
		        $(this).toggleClass('row_selected');  
		    }
		     
		    lastChecked = this;
		    
		    $("#id_goto").val(oTable.fnSettings()._iDisplayStart + $(this).index()+1)
		}	    
	});
}

function initializeMultipleSelectionTool(){
	//Hover on multiple selection tool	
	hiConfig = {
        sensitivity: 3, // number = sensitivity threshold (must be 1 or higher)
        interval: 300, // number = milliseconds for onMouseOver polling interval
        timeout: 800, // number = milliseconds delay before onMouseOut
        over: function(e) {
        	if ($(this).hasClass("row_selected")){
    			$("#multipleSelectionTool").css('left',e.pageX)
    			$("#multipleSelectionTool").css('top',e.pageY)
    			
    			var iRow = $(this).index()
    			$("#multipleSelectionTool").data('row_id', iRow)
    			
    			$("#multipleSelectionTool").fadeIn('slow')
    		}
        }, // function = onMouseOver callback (REQUIRED)
        out: function() { 
        	if ($(this).hasClass("row_selected")){
    			$("#multipleSelectionTool").fadeOut('slow')
    		}
        } // function = onMouseOut callback (REQUIRED)
    }	
	$("tr").hoverIntent(hiConfig)
	
}

function showTableConfig(){
//	Display Div 
	$("#configurationContainer").slideDown('slow')

//	Initialize checkbox in table confirguration container (checked and disable attributes)
	for (var label in jsonTableLayoutConfiguration.columnsLayout){ 
		columnLayoutProperties=jsonTableLayoutConfiguration.columnsLayout[label].columnLayoutProperties
		
//		VISIBLE 
		if (columnLayoutProperties.visible){$("#"+label+"_visible").prop('checked', true)}
		else{$("#"+label+"_visible").prop('checked', false)}
		if (columnLayoutProperties.allowSetVisible){$("#"+label+"_visible").removeProp('disabled')}
		else{$("#"+label+"_visible").prop('disabled', 'disabled')}
		
//		RENDERABLE
		if (columnLayoutProperties.renderable){$("#"+label+"_renderable").prop('checked', true)}
		else{$("#"+label+"_renderable").prop('checked', false)}
		if (columnLayoutProperties.allowSetRenderable){$("#"+label+"_renderable").removeProp('disabled')}
		else{$("#"+label+"_renderable").prop('disabled', 'disabled')}
		
//		EDITABLE
		if (columnLayoutProperties.editable){$("#"+label+"_editable").prop('checked', true)}
		else{$("#"+label+"_editable").prop('checked', false)}
		if (columnLayoutProperties.allowSetEditable){$("#"+label+"_editable").removeProp('disabled')}
		else{$("#"+label+"_editable").prop('disabled', 'disabled')}
	}	
 
}	

function saveTableConfiguration(){
	for (var label in jsonTableLayoutConfiguration.columnsLayout){ 
		columnLayoutProperties=jsonTableLayoutConfiguration.columnsLayout[label].columnLayoutProperties
		columnLayoutProperties.visible=$("#"+label+"_visible").prop('checked')
		//From the image element we get the column header index
		columnId=oTable.fnGetColumnIndex(label)
		oTable.fnSetColumnVis(columnId, columnLayoutProperties.visible);
		
		newValue=$("#"+label+"_renderable").prop('checked')
		if (newValue != columnLayoutProperties.renderable){
			columnLayoutProperties.renderable=newValue
			$("#"+label+"_renderable_icon").find("i").toggleClass("fa-eye").toggleClass("fa-eye-slash")
		}
		
		newValue=$("#"+label+"_editable").prop('checked')
		if (newValue != columnLayoutProperties.editable){
			columnLayoutProperties.editable=newValue
			$("#"+label+"_editable_icon").find("#banElement").toggleClass("fa-times")
		}
		
	}
	
	oTable.fnDraw();
	
	$("#configurationContainer").hide()
}

function multipleEnableDisableImage(mode){
	columnId=oTable.fnGetColumnIndex("enabled")
	columnIdReal=oTable.fnColumnIndexToVisible(columnId)
	
	booleanValue=(mode=='enable')
	integerValue= (booleanValue)? 1 : 0
	
	$(".row_selected").each(function(){
		checkbox_element=$('td:eq('+columnIdReal+')',this).find(":checkbox")
		checkbox_element.prop('checked', booleanValue);
		changes[checkbox_element.attr("id")]=integerValue
	})

	if (!$("#saveButton").hasClass("buttonGreyHovered")){
		$("#saveButton").toggleClass("buttonGreyHovered")
	}
	console.log(changes)
}

function multipleSelect(mode){
	row_id=parseInt($("#multipleSelectionTool").data('row_id'));
	
	$("#data_table tbody tr").each(function(){
		if (mode=='all' || (mode=='from' && $(this).index()>row_id) || (mode=='to' && $(this).index()<row_id) ){
			$(this).addClass("row_selected")
		}	
	})
}
