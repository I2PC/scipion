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

function createProjectForm(title, msg) {

	var msg = msg +"<input type='text' id='newProjName' class='content'/>";
	var funcName = 'createProject';
	
	warningPopup(title, msg, funcName);
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

function deleteProjectForm(projName, title, dialog) {

//	var title = 'Confirm project deletion'
	var msg = "<td class='content' value='"	+ projName +"'>"
			+ "Project " + projName
			+ ", "
			+ dialog 
			+ "</td>";
			
	var funcName = 'deleteProject';
	
	warningPopup(title, msg, funcName);
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


