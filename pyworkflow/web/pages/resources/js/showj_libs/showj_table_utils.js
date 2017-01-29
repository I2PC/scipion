/*****************************************************************************
 *
 * Authors:    Adrian Quintana (aquintana@cnb.csic.es)
 * 			   Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
 * Methods used in the showj_table template 
 * 
 * METHODS LIST:
 * 
 * function renderElements(nRow, aData)
 * -> Render elements on table (Image and checkbox). This method is called from datatable script
 * 	  + nRow is the row to be displayed
 *    + aData the row content
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

/* CONSTANTS */
var COL_RENDER_NONE = 0;
var COL_RENDER_ID = 1;
var COL_RENDER_TEXT = 2;
var COL_RENDER_IMAGE = 3;
var COL_RENDER_CHECKBOX = 4;
var COL_RENDER_SLICE = 5;

/** METHODS ***************************************************************** */



function getColumnsDefinition() {
	var jsonColumnsLayout = jsonTableLayoutConfiguration.columnsLayout;
	var columnId = 0;
	var dataForTable = [];
	
	for ( var i in jsonColumnsLayout) {
		var dataRowForTable = [];
		var columnLayoutConfiguration = jsonColumnsLayout[i];
		
		dataRowForTable["bVisible"] = columnLayoutConfiguration.visible;             
		dataRowForTable["sTitle"] = columnLayoutConfiguration.columnLabel;
		dataRowForTable["sSubTitle"] = columnLayoutConfiguration.columnName;
		dataRowForTable["aTargets"] = [columnId];
//		if (columnLayoutConfiguration.columnLabel == 'weight')
//		{
//			console.log("order by weight")
//			dataRowForTable["aaSorting"] = [columnId, 'asc'];
//		}
		if (columnLayoutConfiguration.columnType == COL_RENDER_CHECKBOX) {
//			console.log("id:"+columnIdReal+" checkbox!")
			dataRowForTable["mRender"] = function (data, type, row){ return colRenderCheckbox(data)}
		
		} 
		if(columnLayoutConfiguration.renderable)
		{
			console.log(columnLayoutConfiguration.columnLabel)
			if (columnLayoutConfiguration.columnType == COL_RENDER_IMAGE) {
				dataRowForTable["mRender"] = function (data, type, row){	return colRenderImage(data)}
			
			// columnType = 5
			} else if (columnLayoutConfiguration.columnType == COL_RENDER_SLICE) {
	//			console.log("id:"+columnIdReal+" render slice!")
				dataRowForTable["mRender"] = function (data, type, row){ return colRenderSlice(data)}
			
			// None
			}
		}
		
		dataForTable.push(dataRowForTable);
		columnId++;
	}
//	console.log("columns Definition:")
//	console.log(dataForTable)
	return dataForTable;
}

function colRenderCheckbox(aData){
	var checkbox_element = '<input type=\"checkbox\" onclick=\"valueChange(this);\" '
	var data = aData
			
	if (data == "True" || data == 1 || data == true) {
		checkbox_element += 'checked>'
	} else {
		checkbox_element += '>'
	}
	return checkbox_element
	
}


function colRenderable(aData, renderFunc){
			
			src = '\"' + getSubDomainURL() + '/render_column/?renderFunc=' + renderFunc;
			src += '&image=' + aData + '\"';
			
			var content_html = '<span style="display:none">' 
					+ aData + '</span>'
					+ '<img class=\"tableImages\"'
					+ ' src=' + src
					+ ' data-real_src=' +src
					+ '/>'
		return content_html;
}


function colRenderImage(aData){
	
	html = colRenderable(aData, "get_image")
	return html
}

function colRenderSlice(aData){
	//var volName = aData.split(".")
	//volName = volName[0] + "_tmp.mrc"
	var volName = aData
	
	src = '\"' + getSubDomainURL() + '/render_column/?renderFunc=get_slice' ;
	src += '&image=' + volName + '\"';
	
	var content_html = '<span style="display:none">' 
			+ aData + '</span>'
			+ '<img class=\"tableImages\" '
			+ ' src=' + src
			+ ' data-real_src=' +src
			+ '/>'
			
	return content_html;
}

function getVisibleColumn(n)
{
	count = 0;
	visibles  = 0 
	for (column in jsonTableLayoutConfiguration.columnsLayout)
	{
		if(jsonTableLayoutConfiguration.columnsLayout[column].visible)
		{
			if(visibles == n)
				return count
			visibles ++
		}
		count ++
	}	
	return -1
}



function arrangeColumns(order, sort)
{
	if(jsonTableLayoutConfiguration.colsOrder)
		orderColumns = jsonTableLayoutConfiguration.colsOrder.trim().split(" ")
	if(jsonTableLayoutConfiguration.colsSortby)
	{
		sortby = jsonTableLayoutConfiguration.colsSortby.trim().split(" ")
		sortbyColumn = sortby[0]
		if (sortby.length == 2)
			direction = sortby[1]
		else
			direction = 'asc'
	}
	
	var j = 0
	for (column in jsonTableLayoutConfiguration.columnsLayout)
	{
		columnLabel = jsonTableLayoutConfiguration.columnsLayout[column].columnLabel
		i = 0
		if(jsonTableLayoutConfiguration.colsOrder)
			for(orderColumn in orderColumns)
			{
				if( columnLabel == orderColumns[orderColumn])
				{
					pos = getVisibleColumn(i)
					if(pos != order[j] && pos != -1)
					{
						pos2 = order.indexOf(j)
						aux = order[pos]
	 					order[pos] = order[pos2]
	 					order[pos2] = aux
					}
				}
				i ++
			}
		if(jsonTableLayoutConfiguration.colsSortby && sortbyColumn == columnLabel)
			sort[0] = [order[j], direction]
		j ++
	}
	
}


function initializeColumnHeader() {
	var headerRow = $("#data_table thead tr")
	var cols = oTable.fnSettings().aoColumns;
	
	for(var x = 0, xLen = cols.length; x < xLen; x++) {
		
		var cellIndex = oTable.fnColumnIndexToVisible(x)
		
		if (cellIndex != null) {
			
			// var elm = $('th:eq(' + cellIndex + ')', headerRow).find("div");
			var elm = $('th:eq(' + cellIndex + ')', headerRow);
			
			// To render is necessary the Name, not the Label
			// var textLabel = cols[x].sTitle
			var textName = cols[x].sSubTitle
			elm.attr("id", textName + "___column_header")

			var properties = jsonTableLayoutConfiguration.columnsLayout[textName]
			var headerIcon = getHeaderWithIcon(textName, properties)
			elm.append(headerIcon)
		
		}
	}
}

function getHeaderWithIcon(text, columnLayoutProperties) {
	var iconElements = ''
	 
	/*
	// allowSetVisible 
	if (columnLayoutProperties.visible && columnLayoutProperties.allowSetVisible){
		iconElements+="<span class=\"css_right\">"
		iconElements+="<a class=\"arrowImage\" onclick=\"enableDisableColumn(event,this)\">"
		iconElements+="<i class='fa fa-times'></i></a></span>"	
	}
	*/

	// allowSetEditable
	if (columnLayoutProperties.allowSetEditable) {
		iconElements += "<span class=\"css_right fa-stack\"><a id=\""
				+ text
				+ "_editable_icon\" onclick=\"enableDisableEditableColumn(event,this);\"><i class=\"fa fa-pencil fa-stack-1x"
		iconElements += "\"></i><i id=\"banElement\" class=\"fa fa-stack-1x "
		if (columnLayoutProperties.renderable) {
			iconElements += "fa-times"
		}
		iconElements += "\"></i></a></span>"
	}
	
	// allowSetRenderable
	if (columnLayoutProperties.allowSetRenderable) {
		iconElements += "<span class=\"css_right\"><a id=\""
				+ text
				+ "_renderable_icon\" onclick=\"javascript:enableDisableRenderColumn(event,this);\"><i class=\"fa "
		if (columnLayoutProperties.renderable) {
			iconElements += "fa-eye-slash"
		} else {
			iconElements += "fa-eye"
		}
		iconElements += "\"></i></a></span>"
	}
	return iconElements;
}


function enableDisableColumn(event, element) {
	// From the image element we get the column header index
	var thCell = $(element).closest('th')
	var iCol = $('th').index(thCell)
	var iCol2 = oTable.fnVisibleToColumnIndex(iCol)

	var bVis = oTable.fnSettings().aoColumns[iCol2].bVisible;

	bvis_var = bVis ? false : true;

	oTable.fnSetColumnVis(iCol2, bvis_var);

	// Update table layout configuration model
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].visible = !bVis

	// Update the session variable
//	 var status = "disable";
//	 if (bvis_var){
//		 status = "enable";
//	 }
//	 updateSession(labelColumn, "visible", status)

	// This will avoid column sorting
	event.stopPropagation()

}

function enableDisableRenderColumn(event, element) {
	// Switch button from on to off or viceversa
	$(element).find("i").toggleClass("fa-eye").toggleClass("fa-eye-slash")

	// From the image element we get the column header index
	var thCell = $(element).closest('th')
	var iCol = $('th').index(thCell)

	// We enable/disable render event
	var nTrs = oTable.fnGetNodes();
	$('td:nth-child(' + (iCol + 1) + ')', nTrs).each(
		function() {
			if ($(element).find("i").hasClass("fa-eye-slash")
					&& $(this).find("img").attr('src') == undefined) {
				/* $(this).find("img").attr("src",$(this).find("img").data('real_src')) */
				setImageSize(false)
				initializeImageLoad(false)
			}
			$(this).find("img").toggle()
			$(this).find("span").toggle()
		})

	// Update table layout configuration model
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].renderable = $(element).find("i").hasClass("fa-eye-slash")

	// Update the session variable
	var status = "enable";
	if ($(element).find("i").attr("class") == "fa fa-eye") {
		status = "disable";
	}
	updateSession(labelColumn, "renderable", status)

	// This will avoid column sorting
	event.stopPropagation()

}

function enableDisableEditableColumn(event, element) {
	$(element).find("#banElement").toggleClass("fa-times")

	// From the image element we get the column header index
	var thCell = $(element).closest('th')
	var iCol = $('th').index(thCell)

	// We enable/disable render event
	var nTrs = oTable.fnGetNodes();
	if (!$(element).find("#banElement").hasClass("fa-times")) {
		$('td:nth-child(' + (iCol + 1) + ')', nTrs).addClass('editable')
		setElementsEditable('td:nth-child(' + (iCol + 1) + ')')
		var status = "enable";
	} else {
		$('td:nth-child(' + (iCol + 1) + ')', nTrs).removeClass('editable')
		$('td:nth-child(' + (iCol + 1) + ')', nTrs).unbind('click');
		var status = "disable";
	}
	$('td:nth-child(' + (iCol + 1) + ')', nTrs).each(function() {
		$(this).toggleClass("editable")
	})

	// Update table layout configuration model
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].editable = !$(
			element).find("#banElement").hasClass("fa-times")

	// Update the session variable
	updateSession(labelColumn, "editable", status)

	// This will avoid column sorting
	event.stopPropagation()

}

function setElementsEditable(elements) {
	if (elements == null) {
		var nTrs = oTable.fnGetNodes();
		for ( var label in jsonTableLayoutConfiguration.columnsLayout) {
			columnId = oTable.fnGetColumnIndex(label)
			columnIdReal = oTable.fnColumnIndexToVisible(columnId)
			if (columnIdReal != null) {
				if (jsonTableLayoutConfiguration.columnsLayout[label].editable) {
					setElementsEditable('td:nth-child(' + (columnIdReal + 1)
							+ ')')
				} else {
					$('td:nth-child(' + (columnIdReal + 1) + ')', nTrs).unbind(
							'click');
				}
			}
		}
	} else {
		$(elements, oTable.fnGetNodes()).editable(function(value, settings) {
			// Get position and real column
			var aPos = oTable.fnGetPosition(this)
			columnIdReal = oTable.fnVisibleToColumnIndex(aPos[1])

			// Get label from column
			var label = oTable.fnSettings().aoColumns[columnIdReal].sTitle

			// ID in the rows? Done in the DataTable declaration
			var aoData = oTable.fnSettings().aoData;
			var nTd = aoData[aPos[0]]._anHidden[oTable.fnGetColumnIndex("id")];
			var rowId = $(nTd).text()

			return (value);
		}, {
			"callback" : function(sValue, y) {
				// Update datatable model
				var aPos = oTable.fnGetPosition(this);
				columnIdReal = oTable.fnVisibleToColumnIndex(aPos[1])
				oTable.fnUpdate(sValue, aPos[0], columnIdReal);
			},
			"height" : "14px"
		});
	}
}

function valueChange(element) {
	var elm = $(element)
	var id = elm.attr("id")
	var element_value = "";
		
	if (elm.is("input:checkbox")) {
		element_value = elm.is(":checked")
		
		//Fix to keep the datatable updated
		if (!element_value){
    		elm.prop("checked", false);
		}else{
			elm.prop("checked", true);
		}
		//Fix to keep the datatable updated
		updateListSession(id, "enabled", "table")
		
	} else {
		element_value = elm.val()
	}
}

function initializeTableWidth() {
	var tableWidth = $("section").width()
	if ($("#table_container").width() > tableWidth) {
		$("#table_container").width(tableWidth);
	}
}

function initializeGoToEvent() {
	$('#id_goto').click(function() {
		updateRowSelection();
	});

	$('#id_goto').change(function() {
		updateRowSelection();
	});
}
function updateRowSelection() {
	oTable.fnDisplayRowWithIndex($("#id_goto").val() - 1)
	var nodes = oTable.fnGetNodes();
	selectedElement = nodes[$("#id_goto").val() - 1]
	selectedElement.click();

	$('html, body').animate({
		scrollTop : $(selectedElement).offset().top
	}, 2000);
}

var lastChecked;
function initializeSelectionRowEvent() {
	$('#data_table').on(
			"click",
			"tr",
			function(event) {

				var start = $('#data_table tbody tr').index(this);
				if (start >= 0) {
					if (!lastChecked) {
						lastChecked = this;
					}

					if (event.shiftKey) {
						var end = $('#data_table tbody tr').index(lastChecked);

						for (i = Math.min(start, end); i <= Math.max(start, end); i++) {
							
							var elm = $('#data_table tbody tr').eq(i);
							if (!elm.hasClass('row_selected')) {
								elm.addClass("row_selected");
								
								/* Elements added to the session list */
								updateListSession(elm.attr("id"), "selected", "table")
							}
						}

						// Clear browser text selection mask
						if (window.getSelection) {
							if (window.getSelection().empty) { // Chrome
								window.getSelection().empty();
							} else if (window.getSelection().removeAllRanges) { // Firefox
								window.getSelection().removeAllRanges();
							}
						} else if (document.selection) { // IE?
							document.selection.empty();
						}
					} else {
						if (!event.metaKey && !event.ctrlKey) {
							$("tr.row_selected").each(function() {
								$(this).removeClass('row_selected');
								
								// remove all from selectedList
								updateListSession($(this).attr("id"), "selected", "table")
							
							});
						}

						$(this).toggleClass('row_selected');
						
						// add/remove to selectedList
						updateListSession($(this).attr("id"), "selected", "table")
					}

					lastChecked = this;

					$("#id_goto").val(
							oTable.fnSettings()._iDisplayStart
									+ $(this).index() + 1)
				}
			});
}

function initializeMultipleSelectionTool() {
	// Hover on multiple selection tool
	var hiConfig = {
		sensitivity : 3, // number = sensitivity threshold (must be 1 or higher)
		interval : 300, // number = milliseconds for onMouseOver polling interval
		timeout : 800, // number = milliseconds delay before onMouseOut
		over : function(e) {
			if ($(this).hasClass("row_selected")) {
				$("#multipleSelectionTool").css('left', e.pageX)
				$("#multipleSelectionTool").css('top', e.pageY)

				var iRow = $(this).index()
				$("#multipleSelectionTool").data('row_id', iRow)
				$("#multipleSelectionTool").fadeIn('slow')
			}
		}, // function = onMouseOver callback (REQUIRED)
		out : function() {
//			if ($(this).hasClass("row_selected")) {
//				$("#multipleSelectionTool").fadeOut('slow')
//			}
		} // function = onMouseOut callback (REQUIRED)
	}
	
	$("tr").hoverIntent(hiConfig)

}

function showTableConfig() {
	// Display Div
	$("#configurationContainer").slideDown('slow')

	// Initialize checkbox in table confirguration container (checked and
	// disable attributes)
	for ( var label in jsonTableLayoutConfiguration.columnsLayout) {
		var columnLayoutProperties = jsonTableLayoutConfiguration.columnsLayout[label]

		// VISIBLE
		if (columnLayoutProperties.visible) {
			$("#" + label + "_visible").prop('checked', true)
		} else {
			$("#" + label + "_visible").prop('checked', false)
		}
		if (columnLayoutProperties.allowSetVisible) {
			$("#" + label + "_visible").removeProp('disabled')
		} else {
			$("#" + label + "_visible").prop('disabled', 'disabled')
		}

		// RENDERABLE
		if (columnLayoutProperties.renderable) {
			$("#" + label + "_renderable").prop('checked', true)
		} else {
			$("#" + label + "_renderable").prop('checked', false)
		}
		if (columnLayoutProperties.allowSetRenderable) {
			$("#" + label + "_renderable").removeProp('disabled')
		} else {
			$("#" + label + "_renderable").prop('disabled', 'disabled')
		}

		// EDITABLE
		if (columnLayoutProperties.editable) {
			$("#" + label + "_editable").prop('checked', true)
		} else {
			$("#" + label + "_editable").prop('checked', false)
		}
		if (columnLayoutProperties.allowSetEditable) {
			$("#" + label + "_editable").removeProp('disabled')
		} else {
			$("#" + label + "_editable").prop('disabled', 'disabled')
		}
	}

}

function saveTableConfiguration() {
	for ( var label in jsonTableLayoutConfiguration.columnsLayout) {
		columnLayoutProperties = jsonTableLayoutConfiguration.columnsLayout[label]
		columnLayoutProperties.visible = $("#" + label + "_visible").prop(
				'checked')

		// From the image element we get the column header index
		columnId = oTable.fnGetColumnIndex(columnLayoutProperties.columnLabel)
		if(columnId != -1)
			oTable.fnSetColumnVis(columnId, columnLayoutProperties.visible);

		var newValue = ""

		newValue = $("#" + label + "_renderable").prop('checked')
		if (newValue != columnLayoutProperties.renderable) {
			columnLayoutProperties.renderable = newValue
			$("#" + label + "_renderable_icon").find("i").toggleClass("fa-eye")
					.toggleClass("fa-eye-slash")
		}

		newValue = $("#" + label + "_editable").prop('checked')
		if (newValue != columnLayoutProperties.editable) {
			columnLayoutProperties.editable = newValue
			$("#" + label + "_editable_icon").find("#banElement").toggleClass(
					"fa-times")
		}

	}

	oTable.fnDraw();

	$("#configurationContainer").hide()
}

function multipleEnableDisableImage(mode) {
	var columnId = oTable.fnGetColumnIndex("enabled")
	var columnIdReal = oTable.fnColumnIndexToVisible(columnId)
	var element_value = "";
	
	$(".row_selected").each(function() {
		var elm = $('td:eq(' + columnIdReal + ')', this).find(":checkbox")
		var id = elm.attr("id")
		
		switch(mode){
			case 'enable':
				if(!elm.is(":checked")){
					elm.prop("checked", true);
					// Update the session list
					updateListSession(id, "enabled", "table")
				}
				break;
		
			case 'disable':
				if(elm.is(":checked")){
					elm.prop("checked", false);
					// Update the session list
					updateListSession(id, "enabled", "table")
				}
				break;
		}
	});
}

function multipleSelect(mode) {
	var row_id = parseInt($("#multipleSelectionTool").data('row_id'));

	$("#data_table tbody tr").each(
			function() {
				if(!$(this).hasClass("row_selected")){
					// Update the session list
					updateListSession($(this).attr("id"), "selected", "table")
					
					if (mode == 'all'
							|| (mode == 'from' && $(this).index() > row_id)
							|| (mode == 'to' && $(this).index() < row_id)) {
						$(this).addClass("row_selected")
						
					}
				}
				
			});
}

function updateSession(label, type, status) {
	var uri = "/update_session_table/?label=" + label + "&type=" + type
		+ "&option=" + status
	var URL = getSubDomainURL() + uri
	$.ajax({
		type : "GET",
		url : URL,
		success : function() {
//			 alert(jsonTableLayoutConfiguration.columnsLayout[label].renderable);
		},
		error: function(){
//			alert("error")
		}
	});
}


