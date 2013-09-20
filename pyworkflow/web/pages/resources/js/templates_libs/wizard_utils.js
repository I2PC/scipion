// *** Common Methods *** //

function selectList(elm, mode) {
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
	}

	// hide the last micrograph
	img.hide();
	img_psd.hide();

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
	
	
	if(mode=="raphael"){
		img.load(putImage(uri, "mic", 250, 250), function() {
			// hide the load img
			load_img.hide();
			// show the new micrograph
			img.show();
		});
	}
	else if(mode=="normal"){
		img.load(img.attr("src", uri), function() {
			// hide the load img
			load_img.hide();
			// show the new micrograph
			img.show();
		});
	}	
}

function selectParticle(elm, mode) {
	var row = $("table#list");
	var oldValue = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != oldValue) {
		// unmark the last option
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', 'no-selected');
	}

	// mark the new option
	row.attr('value', oldValue);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');

	// get the img path
	var path_img = elm.attr('value');
	if (path_img == undefined) {
		path_img = elm.val();
	}

	// load and set the image
	var uri = "/get_image/?image=" + path_img + "&dim=250";
	
	if(mode=="raphael"){
		putImage(uri, "particle", 250, 250);
	}
	else if(mode=="normal"){
		$("img#particle").attr("src", uri);
	}
}

function putImage(url, canvas, width, height) {
	$("div#" + canvas).empty();
	var paper = Raphael(document.getElementById(canvas));
	var img = paper.image(url, 0, 0, width, height);
	
//	img.attr({ "clip-rect": "20,20,30,30" });
//	var rec = paper.rect(0, 0, 250, 250);
//	rec.attr("stroke-width", 0);
//	rec.attr("fill", "url(" + url + ")");	
}

function putCircle(radio, canvas, color) {
	$("#" + canvas).empty();
	var paper = Raphael(document.getElementById(canvas));
	var circle = paper.circle(125, 125, radio);
	circle.attr("fill", color);
	circle.attr("opacity", 0.3);
}

// *** Methods Wizard Downsampling *** //

function compositeDownSampling(elm) {
	$.when(selectList(elm,"normal")).then(previewPsd());
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

// *** Methods Wizard CTF *** //

function compositeCTF(elm) {
	$.when(selectList(elm,"normal")).then(previewPsdFreq());
}

function previewPsdFreq() {

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

	// show the new micrograph
	load.load(putImage(uri, "psd_freq", 250, 250), function() {
		// hide the load img
		load.hide();
	});
}

// *** Methods Wizard Particle Mask *** //

function compositeParticle(elm){
	selectParticle(elm,"raphael");
}

// *** Methods Wizard Volume Mask *** //

function compositeVol(elm){
	selectList(elm,"raphael");
}

// *** Methods Wizard Bandpass filter *** //

function compositeBandpass(elm, low, high, decay) {
	$.when(selectParticle(elm, "normal")).then(previewBandpassFilter(low, high, decay));
}

function previewBandpassFilter(low, high, decay) {
	// get the img path
	var path_img = $("tr#" + $("table#list").attr("value")).attr("value");

	// load and set the image
	var uri = "/get_image_bandpass/?image=" + path_img + "&lowFreq=" + low
			+ "&highFreq=" + high + "&decayFreq=" + decay + "&dim=250";

	$("img#imgFiltered").attr("src", uri);
}

// *** Methods Wizard Gaussian filter *** //

function compositeGaussian(elm) {
	$.when(selectParticle(elm,"normal").then(previewSigma()));
}

function previewSigma() {
	// check values
	var sigma = $("#sigma").val();

	// get the img path
	var path_img = $("tr#" + $("table#list").attr("value")).attr("value");

	// load and set the image
	var uri = "/get_image_gaussian/?image=" + path_img + "&sigmaFreq=" + sigma + "&dim=250";

	$("img#imgFiltered").attr("src", uri);
}


