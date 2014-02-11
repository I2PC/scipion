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
 * Methods used in the showj_base template and all the showj pages which extends from it. 
 * 
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * function reloadImages(forceRecallServer)
 * -> Function to reloadImages from the server while reloading waypoints script
 * 	  If forceRecallServer images are requested even when it was same image (src) 
 * 		
 * function initializeArrowEvents()
 * -> Function to initialize arrow events for all arrowImage class element.
 *    Initialize up and down arrow to the closest input text element (class -> menuInputNumber)
 *    
 * function initializeCheckboxEvents()
 * -> Function to initialize checkbox events so it will reload images on change
 * 
 * function initializeComboBoxEvents(defaultZoom)
 * -> Function to initialize combobox events (block, metadata, volume and axis to reslice) to reload the web on change
 * 
 * function initializeColRowModeEvents(mode)
 * -> Function to initialize events for columns and row mode. 
 * 	  It toggles between: 
 * 		Automatic mode => disable cols & rows textfields
 * 		Manual mode => enable cols & rows textfields and initialize value
 * 	 If mode gallery it hides/shows cols & rows subsection 
 * 
 * function hideShowColsRowsSubSectionMenu()
 * -> Function to hide/show col & rows subsection
 * TODO: This can be done with toggle property
 * 
 * function initializeColsRows()
 * -> Function to initialize Columns and Rows textfields value.
 *    Execute when changing to manual mode
 *    
 * function initializeColsRowsEvents()
 * -> Function to initialize Columns and Rows Events  
 * 
 * function updateGalleryLayout(cols, rows, setElement)
 * -> Calculate real cols & rows depending on cols & rows set by the user
 *    and update gallery layout depending on cols & rows.
 *    setElement is the element that triggered the event and it is mandatory
 *    
 * function updateMainContainerDim(colVal)
 * -> Update main container dimension.
 *    colVal is the number of columns in which images are displayed
 *    
 * function initializeZoomEvents()
 * -> Initialize events on zoom change
 *    Calculate dimension (can be set on perc & pixel) & recall the server	
 * 
 * function initializeImageLoad(forceRecall)
 * -> Function to reloadImages from the server
 * 	  If forceRecallServer images are requested even when it was same image (src)
 * 
 * function saveShowjTable(csrf_token)
 * -> Save table change in metadata
 * 
 * function showHideOptionMenu(){
 * -> Show/Hide option menu
 ******************************************************************************/

 /** METHODS ******************************************************************/

function loadTemplateMode(mode){
	alert("changing mode");
	
	if (mode=="table"){
		template="/showj_table/"
	}
	else if (mode=="gallery"){
		template="/showj_gallery/"
	}
	
	$(function() {
		$.ajax({
			url : template,
			success : function(data) {
					$('div#content_view').html(data);
				}
		});		
	});
}


function reloadImages(forceRecallServer){
	$.waypoints('destroy')
	initializeImageLoad(forceRecallServer)
}

function initializeArrowEvents(){
	$(".arrowImage").click(function(){
		if ($(this).closest(".sectionMenu, .subSectionMenu").find(".menuInputNumber").is("[readonly]")==false){
			if ($(this).attr("src").indexOf("Up") != -1){
				$(this).closest(".sectionMenu, .subSectionMenu").find(".menuInputNumber").val(+$(this).closest(".sectionMenu, .subSectionMenu").find(".menuInputNumber").val()+1)
			}
			else{
				$(this).closest(".sectionMenu, .subSectionMenu").find(".menuInputNumber").val(+$(this).closest(".sectionMenu, .subSectionMenu").find(".menuInputNumber").val()-1)
			}
			$(this).closest(".sectionMenu, .subSectionMenu").find(".menuInputNumber").change();
		}	
	});	
}

function initializeCheckboxEvents(){
	$("#id_mirrorY, #id_applyTransformMatrix, #id_onlyShifts, #id_wrap").click(function(){
		reloadImages(true)
	})
}

function initializeComboBoxEvents(defaultZoom){
	$("#id_blockComboBox, #id_labelsToRenderComboBox, #id_volumesToRenderComboBox, #id_resliceComboBox").change(function(){
		$('#id_zoom').val(defaultZoom)
		$('form#showjForm').submit();
	})
}

function initializeColRowModeEvents(mode){
	$("#colRowMode").click(function(e){
		if ($("#colRowModeImage").attr("src").indexOf("On") != -1){
			/* Change image */
			$("#colRowModeImage").attr("src", $("#colRowModeImage").attr("src").replace("On","Off"))
			
			/* Update hidden variable  */
			$("#id_colRowMode").val("Off");
			
			/* Set textfield to readonly */
			$("#id_cols").attr('readonly', 'True');
			$("#id_rows").attr('readonly', 'True');
			
			/* Allow section div to adjust automatically*/
			$("section").css("width", "");

			/* Display rows & cols menu */
			$("#colsSubSectionMenu, #rowsSubSectionMenu").hide()
			
			e.preventDefault();
		}
		else{
			/* Change image */
			$("#colRowModeImage").attr("src", $("#colRowModeImage").attr("src").replace("Off","On"))
			
			/* Set textfield to editable */
			$("#id_cols").removeAttr('readonly');
			$("#id_rows").removeAttr('readonly');
			
			/* Update hidden variable  */
			$("#id_colRowMode").val("On");
			
			
			if ($("#id_cols").val()==""){
				initializeColsRows()
			}
			else{
				/* Set section div width automatically */
				updateMainContainerDim($("#id_cols").val())
			}gallery
			
			/* Display rows & cols menu */
			$("#colsSubSectionMenu, #rowsSubSectionMenu").show()
			
			e.preventDefault();
		}
	});
	
	if (mode == 'gallery'){
		hideShowColsRowsSubSectionMenu()
	}
}

function hideShowColsRowsSubSectionMenu(){
	if ($("#colRowModeImage").attr("src").indexOf("On") != -1){
		$("#colsSubSectionMenu, #rowsSubSectionMenu").show()
	}
	else{
		$("#colsSubSectionMenu, #rowsSubSectionMenu").hide()
	}
}

function initializeColsRows(){
	colVal = $("#id_cols").val()
	rowVal = $("#id_rows").val()
	if (colVal == "" && rowVal == ""){
		imgContainerWidth =$(".img_container").width() + (parseInt($(".img_container").css("margin-right")) * 2)
		colVal = Math.floor($("section").width() / (imgContainerWidth + 5))
		$("#id_cols").val(colVal)
		rowVal = 0
	}
	
	
	updateGalleryLayout(colVal,rowVal,(colVal != "")?"id_cols":"id_rows")
}

function initializeColsRowsEvents(){
	$("#id_cols, #id_rows").on('click change keyup',function(){
		if (!isNaturalNumber($("#id_cols").val())){$("#id_cols").val(1)}
		if (!isNaturalNumber($("#id_rows").val())){$("#id_rows").val(1)}
		updateGalleryLayout($("#id_cols").val(),$("#id_rows").val(), $(this).attr("id"))
	});
}

function updateGalleryLayout(cols, rows, setElement){
	if (setElement == "id_cols"){
		$("#id_rows").val(Math.ceil($(".img_container").length/cols))
	}
	else{
		$("#id_cols").val(Math.ceil($(".img_container").length/rows))
	}
	updateMainContainerDim($("#id_cols").val())
}

function updateMainContainerDim(colVal){
	imgContainerWidth =$(".img_container").width() + (parseInt($(".img_container").css("padding-right")) * 2)
	sectionWidth = colVal * imgContainerWidth
	$("section").width(sectionWidth)
	
	reloadImages(false)
	
}

function initializeZoomEvents(){
	$('#id_zoom').click(function(){
		/* updateImageDimByUrl(); */
		setImageSize(false)
	});
	
	$('#id_zoom').change(function(){
		/* updateImageDimByUrl(); */
		setImageSize(false)
	});
		
	$('#id_zoom').keyup(function(e){
		var code = e.keyCode || e.which; 
		  	if (code  == 13) {               
				 /* updateImageDimByUrl(); */
				 setImageSize(false)
		  	}
	});	 
}

function initializeImageLoad(forceRecall){
	forceRecallServer = forceRecall
	$('.tableImages').waypoint(function(direction){
		element=$(this)
		if (element.data('real_src') != element.attr("src") || forceRecallServer){ 
			element.attr(
				"src",
				element.data('real_src').concat($("#id_mirrorY").is(':checked')?"&mirrorY":"").
				concat($("#id_applyTransformMatrix").is(':checked')?"&applyTransformMatrix":"").
				concat($("#id_onlyShifts").is(':checked')?"&onlyShifts":"").
				concat($("#id_wrap").is(':checked')?"&wrap":"")
				)
		}
	}, {
		offset: '150%'
	});
}

function saveShowjTable(csrf_token){
	if ($("#saveButton").hasClass("buttonGreyHovered")){
		alert("witol")
		$.ajax({
			type : "POST",
			url : "/save_showj_table/",
			dataType : "json",
			data : {
				changes: JSON.stringify(changes),
				csrfmiddlewaretoken: csrf_token
			}, 
			success : function(json) {
				message = json.message;
				if (message == 'Ok'){
					changes={}
					$("#saveButton").toggleClass("buttonGreyHovered")
				}
				else{
					alert(message)
				}
			}
		});
	}
}

function showHideOptionMenu(){
	$("#optionsButton").toggleClass("buttonGreyHovered")
	$("#optionsList").toggleClass("optionsListOn")
}
