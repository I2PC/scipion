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
 * Methods used in the project template
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * 
 * function createProjectForm()
 * 	->	Dialog web form based in messi.js to verify the option to create a project.
 * 		A name for the project is asked. 
 * 
 * function createProject(elm)
 *  ->	Method to execute the creation for a project.
 *  
 * function deleteProjectForm(projName)
 *  ->	Dialog web form based in messi.js to verify the option to delete a project.
 *  
 * function deleteProject(elm)
 *  ->	Method to execute a delete for a project.
 * 
 ******************************************************************************/

 /** METHODS ******************************************************************/

/*
 * PROJECTS
 */

function createProjectForm(title, msg) {
	var msg = msg +"<input type='text' id='newProjName' class='content'/>";
	var funcName = 'createProject';
	warningPopup(title, msg, funcName);
}

function createProject(elm) {
	var projName = elm.val();
	var URL = getSubDomainURL() + "/create_project/?projectName=" + projName
	$.ajax({
		type : "GET",
		url : URL,
		success : function() {
			var URL2 = getSubDomainURL() + "/project_content/?projectName="+projName
			window.location.href = URL2;
//			window.location.href = "/projects/";
		}
	});
}

function deleteProjectForm(projName, title, dialog) {
	var title = 'Confirm DELETE project ' + projName 
	var msg = "<td class='content' value='"	+ projName +"'>"
//			+ "Project <strong>" + projName + "</strong>"
			+ dialog 
			+ "</td>";
			
	var funcName = 'deleteProject';
	
	warningPopup(title, msg, funcName);
}

function deleteProject(elm) {
	var projName = elm.attr('value');
	var URL = getSubDomainURL() + "/delete_project/?projectName=" + projName
	$.ajax({
		type : "GET",
		url : URL,
		success : function() {
			window.location.href = "/projects/";
		}
	});
}

/* 
 * SERVICE PROJECT
 */ 

function serviceProjForm(){
	var title = 'Confirm project creation'
	var dialog = "<p>Your <strong>Project</strong> will be created.</p>" +
		"<p>This process generates an unique <strong>identification code</strong><p>" +
		"</p>that you have to save to access to the content in other moment.</p><br />" + 
		"<p>If you agree, confirm to generate it.</p>"
		
	var funcName = 'createServProject';
	
	accessPopup(title, dialog, funcName, 'Confirm', 'Cancel');
}

function createServProject(elm) {
	projName = randomString(32, '#aA')
	var URL = getSubDomainURL() + "/create_project/?projectName=" + projName
	$.ajax({
		type : "GET",
		url : URL,
		async: false,
		success : function() {
			var title = "ACCESS CODE"
			var msg = "<p>Your <strong>access code</strong> to the <strong>Project</strong> generated is the next:</p>" +
					"<br /><p><h2>"+ projName  +"</h2></p><br /><p>Please <strong>save this code securely</strong> " +
							"to access to the project in the future.</p>";
			var msg = msg +"<input type='hidden' class='content' value='"+ projName +"'/>";
			var funcName = "goToProject"
				
			accessPopup(title, msg, funcName, 'Go to the project', 'Exit');
		}
	});
}

function goToProjectForm() {
	var title = 'Confirm access code'
	var dialog = "<p>Give me the <strong>identification code</strong> " +
			"to access to your <strong>Project</strong>.</p>" +
			"<p>This code was given to you after the creation of your project.</p>" +
			"<p>If It was forgotten, please contact with the <a style='color:firebrick;' href='#'>Scipion group</a>.</p><br />"
	var msg = dialog +"<p><input type='text' id='code' class='content' style='width:100%;text-align:center;'/></p>";
	var funcName = 'goToProject';
	accessPopup(title, msg, funcName, 'Confirm', 'Cancel');
}

function goToProject(elm) {
	var code = elm.val();
	var URL = getSubDomainURL() + "/check_project_id/?code=" + code
	$.ajax({
		type : "GET",
		url : URL,
		success : function(result) {
			if(result == 1){
				var URL2 = getSubDomainURL() + "/project_content/?projectName="+code+"&mode=service"
				window.location.href = URL2;
			} else {
				var title = "Bad Access"
				var msg = "<p>Wrong <strong>access code</strong> given!</p>"+
					"<p>If It was forgotten, please contact with the <a style='color:firebrick;' href='#'>Scipion group</a>.</p>"
				errorPopup(title, msg)
			}
		}
	});
}


