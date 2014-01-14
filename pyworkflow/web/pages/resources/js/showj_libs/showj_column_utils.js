//Render elements on table (Image and checkbox)
function renderElements(nRow, aData){
	 console.log("row") 
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
		   		if (columnLayoutConfiguration.typeOfColumn == 'image'){
			   		 if (columnLayoutConfiguration.columnLayoutProperties.renderable){
			   			 $('td:eq('+columnIdReal+')', nRow).html( '<span style="display:none">'+aData[columnId]+'</span><img class=\"tableImages\" id=\"'+i+'___'+aData[0]+'\" src=\"/get_image/?image='+aData[columnId]+'\"/>' );
			   		 }
			   		 else{
			   			$('td:eq('+columnIdReal+')', nRow).html( '<span>'+aData[columnId]+'</span><img style="display:none" class=\"tableImages\" id=\"'+i+'___'+aData[0]+'\" src=\"/get_image/?image='+aData[columnId]+'\"/>' );
			   		 }
		   		 }
		   		 else if (columnLayoutConfiguration.typeOfColumn == 'checkbox'){
		   			 checkbox_element = '<input type=\"checkbox\" onclick=\"valueChange(this);\" id=\"'+i+'___'+aData[0]+'\" '
		    			 if (aData[columnId] == 1){checkbox_element += 'checked>'}
		    			 else{checkbox_element += '>'}
		   			 $('td:eq('+columnIdReal+')', nRow).html(checkbox_element);
		   		 }
	   		}	
		 }
		 else{
			 invisibleColumns++;
		 }
		 columnId++;
	 }
}

//Initialize icons in column header (render image, disable column and make column editable)
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

//Add renderable, editable and disable icon to each column
function getHeaderWithIcon(text, columnLayoutProperties){
	
	var iconElements = '' 
	
	if (columnLayoutProperties.allowSetEditable){
		iconElements+="<span class=\"css_right\"><div class=\"iconHeader "
		if (columnLayoutProperties.editable){iconElements+="editTableOn"}else{iconElements+="editTableOff"}
		iconElements+="\" id=\""+text+"_editable_icon\" onclick=\"enableDisableEditableRow(event,this)\"></div></span>"
	}
	
	if (columnLayoutProperties.allowSetRenderable){
		iconElements+="<span class=\"css_right\"><div class=\"iconHeader "
			if (columnLayoutProperties.renderable){iconElements+="renderTableOn"}else{iconElements+="renderTableOff"}
		iconElements+="\" id=\""+text+"_renderable_icon\" onclick=\"enableDisableRenderColumn(event,this)\"></div></span>"
	}
	
	return iconElements; 
}

//Create column definition (class, title, visible) from table layout configuration
function getColumnsDefinition(){
	jsonColumnsLayout=jsonTableLayoutConfiguration.columnsLayout
	columnId=0;
	var dataForTable = []
	for (var i in jsonColumnsLayout){ 
		
		var dataRowForTable = []
		
 		var columnLayoutConfiguration = jsonColumnsLayout[i]
//		if (columnLayoutConfiguration.typeOfColumn == "text" && columnLayoutConfiguration.columnLayoutProperties.editable){
//			dataRowForTable["sClass"]= "editable" 
//		}
//		else if(columnLayoutConfiguration.typeOfColumn == "image"){
//			sClass_Tmp = "imageColumn"
//			if (columnLayoutConfiguration.columnLayoutProperties.renderable){
//				dataRowForTable["sClass"]= sClass_Tmp +" editable"
//			} 
//		} 
		
 		if(!columnLayoutConfiguration.columnLayoutProperties.visible){
			dataRowForTable["bVisible"] = false
		} 
		
		dataRowForTable["sTitle"]= i
		dataRowForTable["aTargets"]=[columnId]
		dataRowForTable["bSortable"]= false
		dataForTable.push(dataRowForTable)
		columnId++;
		
	}
	
	return dataForTable;
}


/* Enable/disable column 
 * Disable datatable column when event triggered from icon disable column header 
 */
function enableDisableColumn(event, element){
	//From the image element we get the column header index
	var thCell= $(element).closest('th')
	var iCol = $('th').index(thCell)
	var iCol2 = oTable.fnVisibleToColumnIndex(iCol)
	
	var bVis = oTable.fnSettings().aoColumns[iCol2].bVisible;
	oTable.fnSetColumnVis( iCol2, bVis ? false : true );

	/*  Update table layout configuration model*/
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].columnLayoutProperties.visible=!bVis
	
	//This will avoid column sorting
	event.stopPropagation() 
}

/* Enable/disable render to column
 * When column cell is image, it will change from img element to path text and viceversa 
 */
function enableDisableRenderColumn(event, element){
	/*  Switch button from on to off or viceversa*/
	$(element).toggleClass("renderTableOn").toggleClass("renderTableOff")
	
	//From the image element we get the column header index
	var thCell= $(element).closest('th')
	var iCol = $('th').index(thCell)
	
	//We enable/disable render event
	var nTrs = oTable.fnGetNodes();
	$('td:nth-child('+(iCol+1)+')', nTrs).each(function(){
		$(this).find("img").toggle()
		$(this).find("span").toggle()
	})

	/*  Update table layout configuration model*/
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].columnLayoutProperties.renderable=$(element).hasClass("renderTableOn")
	
	//This will avoid column sorting
	event.stopPropagation() 
}

/* Enable/disable editable to column
 * When column cell is text or a number, it will change from a text element to an input textfield and viceversa 
 */
function enableDisableEditableRow(event, element){
	$(element).toggleClass("editTableOn").toggleClass("editTableOff")
	
	//From the image element we get the column header index
	var tdCell= $(element).closest('tr').find('td')
	var thCell= $(element).closest('tr').find('th')
		
	//We enable/disable render event
	if ($(element).hasClass("editTableOn")){
		tdCell.addClass('editable')
		setElementsEditable(tdCell)
	}
	else{
		tdCell.removeClass('editable')
		tdCell.unbind('click');
	}
	/* tdCell.each(function(){
		$(this).toggleClass("editable")
	}) */
	
	/*  Update table layout configuration model*/
	var labelColumn = thCell.attr("id").split("___")[0]
	jsonTableLayoutConfiguration.columnsLayout[labelColumn].columnLayoutProperties.editable=$(element).hasClass("editTableOn")
	
	//This will avoid column sorting
	event.stopPropagation()
	
}

/*  Make all the elements editable. Onclick elements will change into textfield (Plugin jeditable)*/
function setElementsEditable(elements){
	/* If no elements are given, all dataset is configured */
	if (elements==null){
		var nTrs = oTable.fnGetNodes();
		for (var label in jsonTableLayoutConfiguration.columnsLayout){ 
			if (jsonTableLayoutConfiguration.columnsLayout[label].columnLayoutProperties.editable){
				setElementsEditable($("#data_table_transpose").find('#'+label+"___column_header").closest('tr').find('td'))
			}
			else{
				$("#data_table_transpose").find('#'+label+"___column_header").closest('tr').find('td').unbind('click');
			}
		}
	}
	else{
		elements.editable(
			/* '../examples_support/editable_ajax.php', */
			function(value, settings){ 
				/* Get label from column*/
				var label = $(this).closest('tr').find('th').attr("id").split("___")[0]
				
				/*NAPA DE LUX TENEMOS QUE PONERLE UN ID A LAS FILAS */
				var aoData = oTable.fnSettings().aoData;
				var nTd = aoData[0]._anHidden[ oTable.fnGetColumnIndex("id") ];
				var rowId =$(nTd).text()
				
				/*Keep changes in global variable */
				if (!$("#saveButton").hasClass("buttonGreyHovered")){
					$("#saveButton").toggleClass("buttonGreyHovered")
				}
				
				/*Keep changes in global variable */
				changes[label+"___"+rowId]=value
			    return(value);
			  },
			{
			 "callback": function( sValue, y ) {
				 
				 /* Get label from column*/
				 var label = $(this).closest('tr').find('th').attr("id").split("___")[0]
				 columnId=oTable.fnGetColumnIndex(label) 
				
				 oTable.fnUpdate( sValue, 0, columnId );
			},
			"height": "14px"
		} );
	}	
}

function valueChange(element){
	element_value = ""
	if ($(element).is("input:checkbox")){
		element_value = $(element).is(":checked")
	}
	else{
		element_value = $(element).val()
	}
	
	/*Keep changes in global variable */
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
	var nodes = oTable.fnGetNodes();
	nodes[$("#id_goto").val()-1].click();
}

/* Control row event selection allowing shift, ctrol and single click */
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
		    	/* $(lastChecked).removeClass('row_selected'); */
		    	if (!event.metaKey && !event.ctrlKey){
		    		$('tr').each(function (){
						$(this).removeClass('row_selected');
					});
			    }
	
		        $(this).toggleClass('row_selected');  
		    }
		     
		    lastChecked = this;
		    
		    $("#id_goto").val($(this).index()+1)
		}	    
	});
}

function initializeMultipleSelectionTool(){
/* 	alert("porakientra")
	$('#data_table').on("hover","tr", function(event) {
		console.log("toloco");
	})
	
	$('td').hover(function(event) {
		alert("td takarrasWithout filter")
		
	},function(event) {
		alert("td takarrasWithout filter")
		
	});*/
	
	$('tr').hover(function(e) {
		if ($(this).hasClass("row_selected")){
			var x= e.pageX
			var y= e.pageY
			
			$("#multipleSelectionTool").css('left',x)
			$("#multipleSelectionTool").css('top',y)
			
			var iRow = $(this).index()
			$("#multipleSelectionTool").data('row_id', iRow)
			
			$("#multipleSelectionTool").fadeIn('slow')
		}	
		
	},function(event) {
		if ($(this).hasClass("row_selected")){
			$("#multipleSelectionTool").delay(2000).fadeOut('slow')
		}
		
	}); 
}

/*Display and initialize div for configuring table layout*/
function showTableConfig(){
	/* Display Div */
	$("#configurationContainer").slideDown()

	/*  Initialize checkbox in table confirguration container (checked and disable attributes)*/
	for (var label in jsonTableLayoutConfiguration.columnsLayout){ 
		columnLayoutProperties=jsonTableLayoutConfiguration.columnsLayout[label].columnLayoutProperties
		
		/* VISIBLE */
		if (columnLayoutProperties.visible){$("#"+label+"_visible").prop('checked', true)}
		else{$("#"+label+"_visible").prop('checked', false)}
		if (columnLayoutProperties.allowSetVisible){$("#"+label+"_visible").removeProp('disabled')}
		else{$("#"+label+"_visible").prop('disabled', 'disabled')}
		
		/* RENDERABLE */
		if (columnLayoutProperties.renderable){$("#"+label+"_renderable").prop('checked', true)}
		else{$("#"+label+"_renderable").prop('checked', false)}
		if (columnLayoutProperties.allowSetRenderable){$("#"+label+"_renderable").removeProp('disabled')}
		else{$("#"+label+"_renderable").prop('disabled', 'disabled')}
		
		/* EDITABLE */
		if (columnLayoutProperties.editable){$("#"+label+"_editable").prop('checked', true)}
		else{$("#"+label+"_editable").prop('checked', false)}
		if (columnLayoutProperties.allowSetEditable){$("#"+label+"_editable").removeProp('disabled')}
		else{$("#"+label+"_editable").prop('disabled', 'disabled')}
	}	
 
}	

/* Save new table configuration and redraw data table */
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
			$("#"+label+"_renderable_icon").toggleClass("renderTableOn").toggleClass("renderTableOff")
		}
		newValue=$("#"+label+"_editable").prop('checked')
		if (newValue != columnLayoutProperties.editable){
			columnLayoutProperties.editable=newValue
			$("#"+label+"_editable_icon").toggleClass("editTableOn").toggleClass("editTableOff")
		}
		
	}
	
	oTable.fnDraw();
	$("#configurationContainer").hide()
}