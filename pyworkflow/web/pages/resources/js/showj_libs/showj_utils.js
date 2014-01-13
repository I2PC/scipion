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

//Initialize block, metadata, volume and axis to reslice combobox to reload the web when changing
function initializeComboBoxEvents(defaultZoom){
	$("#id_blockComboBox, #id_labelsToRenderComboBox, #id_volumesToRenderComboBox, #id_resliceComboBox").change(function(){
		$('#id_zoom').val(defaultZoom)
		$('form#showjForm').submit();
	})
}

//Initialize events for columns and row mode. 
//Automatic mode => disable cols & rows textfields
//Manual mode => enable cols & rows textfields and initialize value
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
			}
			
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

//Initialize Columns and Rows textfields. Execute when changing to manual mode
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

function isNaturalNumber(n) {
    n = n.toString(); // force the value incase it is not
    var n1 = Math.abs(n),
        n2 = parseInt(n, 10);
    
    return !isNaN(n1) && n2 === n1 && n1.toString() === n && n2>0;
}

//Change gallery layout when cols and rows textfield change
function initializeColsRowsEvents(){
	$("#id_cols, #id_rows").on('click change keyup',function(){
		if (!isNaturalNumber($("#id_cols").val())){$("#id_cols").val(1)}
		if (!isNaturalNumber($("#id_rows").val())){$("#id_rows").val(1)}
		updateGalleryLayout($("#id_cols").val(),$("#id_rows").val(), $(this).attr("id"))
	});
}

function updateGalleryLayout(cols, rows, setElement){
	recalculateColsRows(cols,rows, setElement)
	updateMainContainerDim($("#id_cols").val())
}

function recalculateColsRows(cols, rows, setElement){
	if (setElement == "id_cols"){
		$("#id_rows").val(Math.ceil($(".img_container").length/cols))
	}
	else{
		$("#id_cols").val(Math.ceil($(".img_container").length/rows))
	}
}

function updateMainContainerDim(colVal){
	imgContainerWidth =$(".img_container").width() + (parseInt($(".img_container").css("padding-right")) * 2)
	sectionWidth = colVal * imgContainerWidth
	$("section").width(sectionWidth)
	
	reloadImages(false)
	
}

//Initialize events on zoom change (calculate dim (perc & pixel) & recall the server)
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

//Save table change in metadata
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
					/* $("#saveButton").toggleClass("saveOn").toggleClass("saveOff") */
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
