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
 * Methods used in the upload template
 *
 ******************************************************************************/

 /** METHODS ******************************************************************/

// To fill the list of files tab
function doInitFunction(){
	var project_folder = $("#project_folder").val()
	updateListFiles(project_folder)
	var URL = getSubDomainURL() + '/doUpload/'
	$("#uploadForm").submit(function(e) {
		$.ajax({
			url: URL,
			type: 'POST',
			data: new FormData(this),
			processData: false,
			contentType: false,
			dataType: "text",
			async: false,
			success: function(data){
				if(data == "error"){
					errorPopup('Error', "Problem found with the file selected", 0);
				}
				else{
					infoPopup('Success', "The file was uploaded successfully", 0);
				}
			}
		});
		e.preventDefault();
//		$("#progressbar").hide()
		updateListFiles(project_folder)
	});
}

function browseUpload(paramName){
	url_param = "/upload/?mode=service"
	var URL = getSubDomainURL() + url_param
	
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "html",
		async: false,
		success : function(html) {
			new Messi(html, {
				title : 'Select file',
				modal : true,
				buttons : [ {
					id : 0,
					label : 'Upload',
					val : 'Y',
					btnClass : 'fa-cogs',
					btnFunc : 'uploadService',
				}, {
					id : 1,
					label : 'Cancel',
					val : 'C',
					btnClass : 'fa-ban'
				}]
			});
		}
	});
}
 
function launchUpload(){
	var msg = "</td><td class='content' value='"
		msg += "'>The file will be <strong>UPLOADED</strong> into the path <strong>PROJECT FOLDER</strong>. "
		msg += "Do you really want to continue?</td></tr></table>";
	warningPopup('Confirm UPLOAD', msg, 'launchSubmitUpload')
}

function launchSubmitUpload(){
	$('#uploadForm').submit();
}

function uploadService(){
	launchSubmitUpload()
	
	var browser = detectWebBrowser()
	
	if(browser=="chrome" || browser=="ie"){
		var fn = $("input#id_docfile").val().split("fakepath")[1].slice(1).replace(" ", "_")
	}else if("firefox"){
		var fn = $("input#id_docfile").val()
	}else{
		var fn = $("input#id_docfile").val()
		console.log("1:"+fn)
		fn = fn.split("fakepath")[1]
		console.log("2:"+fn)
		fn = fn.slice(1)
		console.log("3:"+fn)
		fn = fn.replace(" ", "_")
		console.log("4:"+fn)
//		var fn = $("input#id_docfile").val().split("fakepath")[1].slice(1).replace(" ", "_")
	}
	
	$("input.upload2").val(fn)

	file = $("#project_folder").val()+ "/Uploads/" + fn
	$("input.upload1").val(file)
}

function updateListFiles(project_folder){
	var project_folder = project_folder + "/Uploads"
	$("tr#listFiles").empty();
	var URL = getSubDomainURL() + "/getPath/?path="+ project_folder
	$.ajax({
		type : "GET",
		url : URL,
		dataType : "json",
		async: false,
		success : function(json) {
			$.each(json, function(key, value) {
				var icon = "<td><img src='"+ value.icon +"' /></td>"
				var url_file = "/get_file/?path="+ project_folder + "/" + value.name +""
				console.log(url_file)
				var name = "<td><a href='" + url_file + "'>"+ value.name +"</a></td>"
				$("tr#listFiles").append("<tr>"+ icon + name + "</tr>")
			});
		}
	});
}