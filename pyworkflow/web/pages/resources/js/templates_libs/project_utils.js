 /**************************************************************************
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
 **************************************************************************/



/**
 * Methods used in the project template
 * 
 * createProjectForm();
 * createProject(elm);
 * deleteProjectForm(projName);
 * deleteProject(elm);
 * 
 **/

function createProjectForm() {

	var html = "Project Name: <input type='text' id='newProjName' class='content'/>";

	new Messi(html, {
		title : 'Enter the project name',
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Ok',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : 'createProject'
		}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'fa-ban'
		} ]
	});
}

function createProject(elm) {
	var projName = elm.val();
	$.ajax({
		type : "GET",
		url : "/create_project/?projectName=" + projName,
		success : function() {
			window.location.href = "/projects/";
		}
	});
}

function deleteProjectForm(projName) {

	var msg = "<table><tr><td><i class=\"fa fa-warning fa-4x\" style=\"color:#fad003;\"></i>"
			+ "</td><td class='content' value='"
			+ projName
			+ "'>Are you sure to <strong>DELETE</strong> project <strong>"
			+ projName
			+ "</strong> and all its <strong>DATA</strong></td></tr></table>";

	new Messi(msg, {
		title : 'Confirm project deletion',
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : 'deleteProject'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'fa-ban'
		} ]
	});
}

function deleteProject(elm) {
	var projName = elm.attr('value');
	$.ajax({
		type : "GET",
		url : "/delete_project/?projectName=" + projName,
		success : function() {
			window.location.href = "/projects/";
		}
	});
}


