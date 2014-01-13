function initRenderFunction(labelToRender){
	renderFunc=jsonTableLayoutConfiguration.columnsLayout[labelToRender].columnLayoutProperties.renderFunc
	extraRenderFunc=jsonTableLayoutConfiguration.columnsLayout[labelToRender].columnLayoutProperties.extraRenderFunc

	console.log(jsonTableLayoutConfiguration.columnsLayout) 
	
	if (renderFunc!=""){ 
		$(".tableImages").each(function(){
		   	var newSrc = $(this).data("real_src").replace(/(renderFunc=).*?(&)/,'$1' + renderFunc + '$2');
		   	if (extraRenderFunc!=""){newSrc = newSrc.concat("&"+extraRenderFunc)}
			$(this).data("real_src", newSrc);
			console.log("newSrc",newSrc)
	   	})
	}
}

function initializeGoToEvents(){
	$('#id_goto').click(function(){
		updateSelection();	
	});
	
	$('#id_goto').change(function(){
		updateSelection();
	});
		
	$('#id_goto').keyup(function(e){
		var code = e.keyCode || e.which; 
		  	if (code  == 13) {               
		  		updateSelection();
		  	}
	});	 
}

function updateSelection(){
	$(".img_container").each(function(){
		$(this).removeClass('image_selected')
	});
	selectedElement=$("#img_container___"+$('#id_goto').val())
	selectedElement.addClass("image_selected")
	
	$('html, body').animate({
        scrollTop: selectedElement.offset().top
    }, 2000);
}

//Show enabled/disabled image and title on hover
function initializeImageEvents(hasEnabledColumn){
	//Enable/Disable image
	$(".img_container").hover(
		function(e){
			if (!$(this).find(".enabledGallery").hasClass("selected")){
				$(this).find(".enabledGallery").fadeIn('slow')
			}
		},
		function(){
			if (!$(this).find(".enabledGallery").hasClass("selected")){
				$(this).find(".enabledGallery").fadeOut('fast')
			}
		}) 
	
		
	if (hasEnabledColumn){	
			
	//Hover on multiple selection tool	
	hiConfig = {
        sensitivity: 3, // number = sensitivity threshold (must be 1 or higher)
        interval: 300, // number = milliseconds for onMouseOver polling interval
        timeout: 800, // number = milliseconds delay before onMouseOut
        over: function(e) {
        	if ($(this).hasClass("image_selected") && $("#multipleSelectionTool").css('display')!="block"){
        		hoverElement=$(this).attr("id")
									
				$("#multipleSelectionTool").css('left',e.pageX)
				$("#multipleSelectionTool").css('top',e.pageY)
				$("#multipleSelectionTool").data('img_container_id',hoverElement)
				$("#multipleSelectionTool").fadeIn('slow')
			}
        }, // function = onMouseOver callback (REQUIRED)
        out: function() { 
        	if ($(this).hasClass("image_selected") && hoverElement==$(this).attr("id")){
				$("#multipleSelectionTool").fadeOut('slow')
			}
        } // function = onMouseOut callback (REQUIRED)
    }	
	$(".img_container").hoverIntent(hiConfig)
	}
	
	//Selection tool
	$(".img_container").click(
		function(e){
			element_id=parseInt($(this).attr("id").split("___")[1]);

			if (e.shiftKey){
				prev_element_id=parseInt($('#id_goto').val())
				
				if (prev_element_id < element_id){
					initialIndex = prev_element_id
					endIndex = element_id
				}
				else{
					initialIndex = element_id
					endIndex = prev_element_id
				}
				
				for (var x=initialIndex; x<=endIndex; x++){
					$("#img_container___"+x).addClass("image_selected")
				}
			}
			else{
				if (!e.ctrlKey && !e.metaKey){
					$(".img_container").each(function(){
						$(this).removeClass('image_selected')
					});	
				}
				$("#img_container___"+element_id).toggleClass("image_selected")
			}
			$('#id_goto').val(element_id);
		})	
}

function initializeWindowLoadEvents(){
	$(window).load(
	function(){
		if ($("#id_colRowMode").val() == 'On' && $("#id_mode").val() == 'gallery'){
			updateMainContainerDim($("#id_cols").val())
		}
	})
}

function enableDisableImage(element, enableDisable){
	if (enableDisable=='enable'){
		$(element).removeClass("selected")
	}
	else if (enableDisable=='disable'){
		$(element).addClass("selected")
	}
	else{
		$(element).toggleClass("selected")
	}
	
	
	if ($(element).hasClass("selected")){element_value = 0}
	else{element_value = 1}
	
//		Keep changes in global variable
	if (!$("#saveButton").hasClass("buttonGreyHovered")){
		$("#saveButton").toggleClass("buttonGreyHovered")
	}
		
	changes[$(element).attr("id")]=element_value
}
function multipleEnableDisableImage(mode){
	if (mode=='enable'){
		$(".image_selected").each(function(){
			enableDisableImage($(this).find(".enabledGallery"), 'enable')
			$(this).find(".enabledGallery").fadeOut('fast')
		})
	}
	else{
		$(".image_selected").each(function(){
			enableDisableImage($(this).find(".enabledGallery"), 'disable')
			$(this).find(".enabledGallery").fadeIn('fast')
		})
	}
}
function multipleSelect(mode){
	element_id=parseInt($("#multipleSelectionTool").data('img_container_id').split("___")[1]);
	$(".img_container").each(function(){
		if (mode=='all' || (mode=='from' && $(this).attr('id').split("___")[1]>element_id) || (mode=='to' && $(this).attr('id').split("___")[1]<element_id) ){
			$(this).addClass("image_selected")
		}	
	})
}