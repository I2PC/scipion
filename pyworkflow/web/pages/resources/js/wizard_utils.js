 /*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
 * Methods to manage wizards.
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * 
 * // *** Common Methods *** //
 * function selectList(elm, mode)
 * 	->	Function to select an image from a list and be loaded with a 
 * 		specific mode (Normal or using the library Raphael.js).
 * 
 * function selectParticle(elm, mode)
 * 	->	Function to select a particle from a set and be loaded with a 
 * 		specific mode (Normal or using the library Raphael.js).
 * 
 * function putImage(url, canvas, width, height)
 * 	->	Overwrite a canvas space with a new image using the library Raphael.js
 * 
 * function putCircle(radio, canvas, color)
 * 	->	Draw a circle over a canvas using the library Raphael.js
 * 
 * function previewPsd()
 * 	->	This function obtain an PSD image from another image. The downsampling
 * 		factor is used here to compute the new image.
 *
 * // *** Methods Wizard Downsampling *** //
 * function previewDownSampling(elm)
 * 	->	Function composite structured in two steps:
 * 			1. Select an image from a list.
 * 			2. To obtain the downsampled image from the image selected.
 * 
 * // *** Methods Wizard CTF *** //
 * function previewCTF(elm)
 * 	->	Function composite structured in two steps:
 * 			1. Select an image from a list.
 * 			2. To obtain a preview image using the method previewPSD() to 
 * 			visualize the CTF.

 * // *** Methods Wizard Particle Mask *** //
 * function compositeParticle(elm)
 * 	->	Function to select a particle from a set and be loaded using the 
 * 		library Raphael.js. Method specific for particles.
 * 
 * // *** Methods Wizard Volume Mask *** //
 * function compositeVol(elm)
 * 	->	Function to select an image from a list and be loaded using the library
 * 		Raphael.js. Method specific for volumes.
 * 
 * // *** Methods Wizard Bandpass filter *** //
 * function compositeBandpass(elm, low, high, decay)
 * 	->	Function composite structured in two steps:
 * 			1. Select a particle from a list.
 * 			2. Apply a bandpass filter to the image based in three parameters and
 * 			show the preview.
 * 
 * function previewBandpassFilter(low, high, decay)
 * 	->	This function get a image using a bandpass filter based in the three 
 * 		parameters: low frequency, high frequency and decay.
 * 
 * // *** Methods Wizard Gaussian filter *** //
 * function compositeGaussian(elm)
 * 	->	Function composite structured in two steps:
 * 			1. Select a particle from a list.
 * 			2. Apply a gaussian filter to the image based in the sigma parameter.
 * 
 * function previewSigma()
 * 	->	This function get a image using a gaussian filter based in the sigma
 * 		parameter.
 * 
 * // *** Methods Wizard Spider Particle filter *** //
 * function compositeSpiderFilter(elm, filterType, filterMode, usePadding)
 * 	->	Function composite structured in two steps:
 * 			1. Select a particle from a list.
 * 			2. Apply a spider filter to the image based.
 * 
 * function previewSpiderFilter(filterType, filterMode, usePadding)
 * 	->	This function get a image using a spider filter with some parameters.
 *
 * function convertBandPass(low, high, decay, samplingRate)
 * 	->	This function convert band pass parameters to A.
 *
 ******************************************************************************/

/** METHODS ******************************************************************/

// *** Common Methods *** //
function selectList(elm, mode) {
	/*
	 * Function to select an image from a list and be loaded with a 
	 * specific mode (Normal or using the library Raphael.js) 
	 */
	var row = $("table#data");
	var oldValue = elm.attr('id');

	var img = $("img#mic");
	var img_psd = $("img#psd");

	var load_img = $("img#loadingMic");

	if (row.attr('value') != undefined && row.attr('value') != oldValue) {
		// unmark the last option
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', '');
		rowOld.attr('class', 'no-selected');
	}

	// hide the last micrograph
	img.hide();
	img_psd.hide();

	// mark the new option
	row.attr('value', oldValue);
	elm.attr('style', 'background-color: #F3CBCB;');
	elm.attr('class', 'selected');

	// set loading img
	load_img.show();

	// get the img path
	var path_img = elm.attr('value');
	if (path_img == undefined) {
		path_img = elm.val();
	}

	// load and set the image
	var URL = getSubDomainURL() + "/get_image/?image=" + path_img + "&dim=250"
	
	if(mode=="raphael"){
		img.load(putImage(URL, "mic", 250, 250), function() {
			// hide the load img
			load_img.hide();
			// show the new micrograph
			img.show();
		});
	}
	else if(mode=="normal"){
		img.load(img.attr("src", URL), function() {
			// hide the load img
			load_img.hide();
			// show the new micrograph
			img.show();
		});
	}	
}

function selectParticle(elm, mode) {
	/*
	 * Function to select a particle from a set and be loaded with a 
	 * specific mode (Normal or using the library Raphael.js)
	 */
	var row = $("table#data");
	var oldValue = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != oldValue) {
		// unmark the last option
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', '');
		rowOld.attr('class', 'no-selected');
	}

	// mark the new option
	row.attr('value', oldValue);
	elm.attr('style', 'background-color: #F3CBCB;');
	elm.attr('class', 'selected');

	// get the img path
	var path_img = elm.attr('value');
	if (path_img == undefined) {
		path_img = elm.val();
	}

	// load and set the image
	var URL = getSubDomainURL() + "/get_image/?image=" + path_img + "&dim=250"
	
	if(mode=="raphael"){
		putImage(URL, "particle", 250, 250);
	}
	else if(mode=="normal"){
		$("img#particle").attr("src", URL);
	}
}

function putImage(url, canvas, width, height, mode) {
	/*
	 * Overwrite a canvas space with a new image using the library Raphael.js 
	 */
	 	 
	$("div#" + canvas).empty();
	var paper = Raphael(document.getElementById(canvas));
	var img = paper.image(url, 0, 0, width, height);
	
}

function putImageSrc(url, canvas) {
	var img = $("img#" + canvas)
	img.attr("src", url)	
}

function putCircle(radio, canvas, color) {
	$("#" + canvas).empty();
	var paper = Raphael(document.getElementById(canvas));
	var circle = paper.circle(125, 125, radio);
	circle.attr("fill", color);
	circle.attr("opacity", 0.3);
}

function previewPsd(){
	/*
	 * This function obtain an PSD image from another image. The downsampling
	 * factor is used here to compute the new image.
	 */
	
	// check downsampling is a number
	var downsampling = $("input#downsampling").val();
	
	if (downsampling == undefined) {
		downsampling = 1.0;
	}

	// get the img path
	var path_img = $("tr#" + $("table#data").attr("value")).attr("value");
	
	// set loading img
	var img_load = $("img#loadingPsd");
	img_load.show();

	// load and set the image
	var uri = "/get_image_psd/?image=" + path_img + "&downsample="
			+ downsampling + "&dim=250"
	var URL = getSubDomainURL() + uri
			
	return [URL, img_load];
}

// *** Methods Wizard Downsampling *** //
function previewDownsampling(elm) {
	/*
	 * Function composite structured in two steps:
	 *	1. Select an image from a list.
	 * 	2. To obtain the downsampled image from the image selected.
	 */
	
	if (elm == undefined){
		elm = $("tr.selected")
	}
	else{
		selectList(elm,"normal")
	}
	
	var img = $("img#psd");
	img.hide();
	
	var res = previewPsd();
	var uri = res[0];
	var img_load = res[1];
	
	img.load(img.attr("src", uri), function() {
		// hide the load img
		img_load.hide();
		// show the new micrograph
		img.show();
	});
	
	
//	$.when(selectList(elm,"normal")).then(
//		function(){
//			var img = $("img#psd");
//			img.hide();
//			
//			var res = previewPsd();
//			var uri = res[0];
//			var img_load = res[1];
//			
//			img.load(img.attr("src", uri), function() {
//				// hide the load img
//				img_load.hide();
//				// show the new micrograph
//				img.show();
//			});
//		}
//	);
		
	 
}

// *** Methods Wizard CTF *** //
function previewCTF(elm) {
	/*
	 * Function composite structured in two steps:
	 * 	1. Select an image from a list.
	 * 	2. To obtain a preview image using the method previewPSD() to 
	 * 	visualize the CTF.
	 */
	$.when(selectList(elm,"normal")).then(
		function(){
			var img = $("img#psd");
			img.hide();
			
//			var uri, img_load = previewPsd()
			var res = previewPsd();
			var uri = res[0];
			var img_load = res[1];
			
			// show the new micrograph
			img.load(putImage(uri, "psd_freq", 250, 250), function() {
				// hide the load img
				img_load.hide();
			});
		}
	);
}

function doCTF() {
	/*
	 * 	To obtain a preview image using the method previewPSD() to 
	 * 	visualize the CTF.
	 */
	var img = $("img#psd");
	img.hide();
	
//			var uri, img_load = previewPsd()
	var res = previewPsd();
	var uri = res[0];
	var img_load = res[1];
	
	// show the new micrograph
	img.load(putImage(uri, "psd_freq", 250, 250), function() {
		// hide the load img
		img_load.hide();
	});
}

// *** Methods Wizard Particle Mask *** //
function compositeParticle(elm){
	/*
	 * Function to select a particle from a set and be loaded using the 
	 * library Raphael.js. Method specific for particles.
	 */
	selectParticle(elm,"raphael");
}

// *** Methods Wizard Volume Mask *** //
function compositeVol(elm){
	/*
	 * Function to select an image from a list and be loaded using the library
	 * Raphael.js. Method specific for volumes.
	 */
	selectList(elm,"raphael");
}

// *** Methods Wizard Bandpass filter *** //
function compositeBandpass(elm, low, high, decay) {
	/*
	 * Function composite structured in two steps:
	 * 	1. Select a particle from a list.
	 * 	2. Apply a bandpass filter to the image based in three parameters and
	 * 	show the preview.
	 */
	$.when(selectParticle(elm, "normal")).then(previewBandpassFilter(low, high, decay));
}

function compositeBandpassVol(elm, low, high, decay) {
	/*
	 * Function composite structured in two steps:
	 * 	1. Select a volume from a list.
	 * 	2. Apply a bandpass filter to the image based in three parameters and
	 * 	show the preview.
	 */
	$.when(selectList(elm, "normal")).then(previewBandpassFilter(low, high, decay));
}

function previewBandpassFilter(low, high, decay) {
	// get the img path
	var path_img = $("tr#" + $("table#data").attr("value")).attr("value");

	// load and set the image
	var uri = "/get_image_bandpass/?image=" + path_img + "&lowFreq=" + low
			+ "&highFreq=" + high + "&decayFreq=" + decay + "&dim=250";
	var URL = getSubDomainURL() + uri
	
	$("img#imgFiltered").attr("src", URL);
}

// *** Methods Wizard Gaussian filter *** //

function compositeGaussian(elm) {
	/*
	 * Function composite structured in two steps:
	 * 	1. Select a particle from a list.
	 * 	2. Apply a gaussian filter to the image based in the sigma parameter.
	 */
	selectParticle(elm,"normal")
	previewSigma()
//	$.when(selectParticle(elm,"normal").then(previewSigma()));
}

function compositeGaussianVol(elm) {
	/*
	 * Function composite structured in two steps:
	 * 	1. Select a volume from a list.
	 * 	2. Apply a gaussian filter to the image based in the sigma parameter.
	 */
	selectList(elm,"normal")
	previewSigma()
//	$.when(selectList(elm,"normal").then(previewSigma()));
}

function previewSigma() {
	/*
	 * This function get a image using a gaussian filter based in the sigma
	 * parameter.
	 */
	// check values
	var sigma = $("#sigma").val();

	// get the img path
	var path_img = $("tr#" + $("table#data").attr("value")).attr("value");

	// load and set the image
	var uri = "/get_image_gaussian/?image=" + path_img + "&sigmaFreq=" + sigma + "&dim=250";
	var URL = getSubDomainURL() + uri
	
	$("img#imgFiltered").attr("src", URL);
}

// *** Methods Wizard Spider Particle filter *** //
function compositeSpiderFilter(elm, filterType, filterMode, usePadding){
	/*
	 * Function composite structured in two steps:
	 * 	1. Select a particle from a list.
	 * 	2. Apply a spider filter to the image based.
	 * 
	 */
	selectParticle(elm,"raphael");
	previewSpiderFilter(filterType, filterMode, usePadding);
}

function previewSpiderFilter(filterType, filterMode, usePadding) {
	/*
	 * This function get a image using a spider filter with some parameters.
	 */
	// get the img path
	var path_img = $("tr#" + $("table#data").attr("value")).attr("value");
	
	// load and set the image
	var uri = "/get_image_filter_spider/?image=" + path_img + "&filterType=" + filterType + "&dim=250"+ "&filterMode="+filterMode+"&usePadding="+usePadding;
	var URL = getSubDomainURL() + uri
	
	// check values
	if (filterType < 2){
		var radius = $('input#radius_val').val();
		URL += "&radius="+radius;
	} else {
		var highFreq = $('input#high_val').val();
		var lowFreq = $('input#low_val').val();
		URL += "&highFreq="+highFreq+"&lowFreq="+lowFreq;
		
		if (filterType == 2){
			var temperature = $('input#temp_val').val();
			URL += "&temperature="+temperature;
		}
	}
	putImage(URL, "filteredParticle", 250, 250);
}


function previewSpiderCustomMask(path, radius1, sdFactor, radius2, maskThreshold) {
	/*
	 * This function get a image using a spider custom mask with some parameters.
	 */
	
	// load and set the image
	var URL = getSubDomainURL() + "/run_custom_mask_spider/?" +
				"image=" + path +
				"&radius1=" + radius1 +
				"&sdFactor=" + sdFactor +
				"&radius2=" + radius2 +
				"&maskThreshold="+ maskThreshold


	 //	Generate the images
	$.ajax({
		type : "GET",
		url : URL,
		async : false
	});

	var canvasList = ['image1', 'image2', 'image3',
	                  'image4', 'image5', 'image6',
	                  'image7', 'image8']

	// Put the images in his canvas associated
	var URL2 = getSubDomainURL() + "/get_image/?image=" + path + "&dim=100";
	putImageSrc(URL2, canvasList[0]);

    // Pass current time as prefix to avoid that image is not refresh when cached by browser
	for (var i=1;i < canvasList.length;i++){
		URL2 = getSubDomainURL() + "/get_image/?image=" + i + "@stkmask.stk&dim=100&prefix=" + new Date().getTime();
		putImageSrc(URL2, canvasList[i]);
	}
}

function convertBandPass(low, high, decay, samplingRate) {
    if (low != 0){
        low = samplingRate/low;
    }
    if (high != 0){
        high = samplingRate/high;
    }
    if (decay != 0){
        decay = samplingRate/decay;
    }

    return [low, high, decay];
}
