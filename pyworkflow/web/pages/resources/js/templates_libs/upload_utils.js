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
 * Methods used in the upload template
 *
 ******************************************************************************/

 /** METHODS ******************************************************************/

// To fill the list of files tab
$(document).ready(function() {
	var project_folder = $("#project_folder").val()
	updateListFiles(project_folder)
	var URL = getSubDomainURL() + '/doUpload/'
	$("#uploadForm").submit(function(e) {
		$.ajax({
			url: URL,
	        type: 'POST',
	        data: new FormData( this ),
	        processData: false,
	        contentType: false,
	        async: false
		});
		e.preventDefault();
		infoPopup('Success', "The file was uploaded successfuly",0);
//		$("#progressbar").hide()
		updateListFiles(project_folder)
	});
});
 
 
function launchUpload(){
	var msg = "</td><td class='content' value='"
		msg += "'>The file will be <strong>UPLOADED</strong> into the path <strong>PROJECT FOLDER</strong>. "
		msg += "Do you really want to continue?</td></tr></table>";
	warningPopup('Confirm UPLOAD',msg, 'launchSubmitUpload')
}

function launchSubmitUpload(){
//	$("#progressbar").show()
	$('#uploadForm').submit();
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