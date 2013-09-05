function selectList(elm) {
	var row = $("table#list");
	var oldValue = elm.attr('id');

	var img = $("img#mic");
	var img_psd = $("img#psd");

	var load_img = $("img#loadingMic");

	if (row.attr('value') != undefined && row.attr('value') != oldValue) {
		// unmark the last option
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', 'no-selected');

		// hide the last micrograph
		img.hide();
		img_psd.hide();
	}

	// mark the new option
	row.attr('value', oldValue);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');

	// set loading img
	load_img.show();

	// get the img path
	var path_img = elm.attr('value');
	if (path_img == undefined) {
		path_img = elm.val();
	}

	// load and set the image
	var uri = "/get_image/?image=" + path_img + "&dim=250";

	img.load(img.attr("src", uri), function() {
		// hide the load img
		load_img.hide();
		// show the new micrograph
		img.show();
	});
}

function previewPsd() {
	var img = $("img#psd");

	img.hide();

	// check downsampling is a number
	var downsampling = $("input#downsampling").val();
	if (downsampling == undefined) {
		downsampling = 1.0;
	}

	// get the img path
	var path_img = $("tr#" + $("table#list").attr("value")).attr("value");

	var load = $("img#loadingPsd");

	// set loading img
	load.show();

	// load and set the image
	var uri = "/get_image_psd/?image=" + path_img + "&downsample="
			+ downsampling + "&dim=250";

	img.load(img.attr("src", uri), function() {
		// hide the load img
		load.hide();
		// show the new micrograph
		img.show();
	});
}

function compositePreview(elm) {
	$.when(selectList(elm)).then(previewPsd());

	//	selectList(elm);
	//	previewPsd();
}

//$("div.highFreq").noUiSlider({
//	range : [ 20, 100 ],
//	start : [ 40, 80 ],
//	step : 20,
//	slide : function() {
//		var values = $(this).val();
//		$("span").text(values[0] + " - " + values[1]);
//	}
//});
