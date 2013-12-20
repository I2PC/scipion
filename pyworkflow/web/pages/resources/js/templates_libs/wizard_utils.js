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
 *  e-mail address 'jmdelarosa@cnb.csic.es'
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
 ******************************************************************************/

/** METHODS ******************************************************************/

// *** Common Methods *** //
function selectList(elm, mode) {
	/*
	 * Function to select an image from a list and be loaded with a 
	 * specific mode (Normal or using the library Raphael.js) 
	 */
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
	/*
	 * Function to select a particle from a set and be loaded with a 
	 * specific mode (Normal or using the library Raphael.js)
	 */
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
	/*
	 * Overwrite a canvas space with a new image using the library Raphael.js 
	 */
	$("div#" + canvas).empty();
	var paper = Raphael(document.getElementById(canvas));
	var img = paper.image(url, 0, 0, width, height);
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
	var path_img = $("tr#" + $("table#list").attr("value")).attr("value");

	// set loading img
	var load = $("img#loadingPsd");
	load.show();

	// load and set the image
	var uri = "/get_image_psd/?image=" + path_img + "&downsample="
			+ downsampling + "&dim=250";
	
	return uri, load
}

// *** Methods Wizard Downsampling *** //
function previewDownSampling(elm) {
	/*
	 * Function composite structured in two steps:
	 *	1. Select an image from a list.
	 * 	2. To obtain the downsampled image from the image selected.
	 */
	$.when(selectList(elm,"normal")).then(
		function(){
			var img = $("img#psd");
			img.hide();
			var uri, load = previewPsd()
			img.load(img.attr("src", uri), function() {
				// hide the load img
				load.hide();
				// show the new micrograph
				img.show();
			});
		}
	);
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
			var uri, load = previewPsd()
			// show the new micrograph
			load.load(putImage(uri, "psd_freq", 250, 250), function() {
				// hide the load img
				load.hide();
			});
		}
	);
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
	/*
	 * Function composite structured in two steps:
	 * 	1. Select a particle from a list.
	 * 	2. Apply a gaussian filter to the image based in the sigma parameter.
	 */
	$.when(selectParticle(elm,"normal").then(previewSigma()));
}

function previewSigma() {
	/*
	 * This function get a image using a gaussian filter based in the sigma
	 * parameter.
	 */
	// check values
	var sigma = $("#sigma").val();

	// get the img path
	var path_img = $("tr#" + $("table#list").attr("value")).attr("value");

	// load and set the image
	var uri = "/get_image_gaussian/?image=" + path_img + "&sigmaFreq=" + sigma + "&dim=250";

	$("img#imgFiltered").attr("src", uri);
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
	var path_img = $("tr#" + $("table#list").attr("value")).attr("value");

	// load and set the image
	var uri = "/get_image_filter_spider/?image=" + path_img + "&filterType=" + filterType + "&dim=250"+ "&filterMode="+filterMode+"&usePadding="+usePadding;

	// check values
	if (filterType < 2){
		var radius = $('input#radius_val').val();
		uri += "&radius="+radius;
	} else {
		var highFreq = $('input#high_val').val();
		var lowFreq = $('input#low_val').val();
		uri += "&highFreq="+highFreq+"&lowFreq="+lowFreq;
		
		if (filterType == 2){
			var temperature = $('input#temp_val').val();
			uri += "&temperature="+temperature;
		}
	}
	putImage(uri, "filteredParticle", 250, 250);
}


