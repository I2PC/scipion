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
			var URL2 = getSubDomainURL() + "/projects/"
			window.location.href = URL2;
		}
	});
}

/* 
 * SERVICE PROJECT
 */ 

function serviceProjForm(){
	var title = 'Project creation'
	var dialog = "<p>Your <strong>Project</strong> will be created.<br />" +
        "This process generates a unique <strong>identification code</strong><br />" +
        "that you will have to save to access the content in the future.</p>" +
        "<p><br /></p>";
	
    dialog += "<p>Confirm to generate it.</p>";

	var funcName = 'createServProject';

	accessPopup(title, dialog, funcName, 'Confirm', 'Cancel');
}

function serviceTestDataForm(){
	var title = 'Test data'
	var dialog = ""
		
	dialog += "<p>You have two options to use the <strong>Test data</strong> :<br />" +
        "1.- <strong>Create a project</strong> with the <strong>Test data</strong> already imported inside.<br />" +	
		"2.- <strong>Download</strong> the files chosen to be imported manually by the user into a project created.</p>";
	dialog += "<br />";
	dialog += '<div id="testData">';
	dialog += "<p>Select the <strong>test data</strong>:</p>";
	dialog += '&nbsp;&nbsp;&nbsp;<input type="radio" name="data" value="groel" checked>';
	dialog += '&nbsp;&nbsp;' + getRefTestData("groel");
	dialog += '<br />';
	dialog += '&nbsp;&nbsp;&nbsp;<input type="radio" name="data" value="bpv">';
	dialog += '&nbsp;&nbsp;' + getRefTestData("bpv");
	dialog += '<br />';
	dialog += '&nbsp;&nbsp;&nbsp;<input type="radio" name="data" value="ribosome">';
	dialog += '&nbsp;&nbsp;' + getRefTestData("ribosome");
	dialog += '<br />';
	dialog += "</div>";
	dialog += "<br />";
    dialog += "<p>After select, choose your option.</p>";

    
    var btn1 = 'Create project'
	var ico1 = 'fa-check'
	var funcName1 = 'createServProject';
		
	var btn2 = 'Download'
	var ico2 = 'fa-download';
	var funcName2 = 'downloadTestdata';
	
	accessPopupOpt(title, dialog, 
					 btn1, ico1, funcName1, 
					 btn2, ico2, funcName2, 
					 "Cancel")
}

function goExampleForm(){
	var title = 'Example projects'
	var dialog = ""

	dialog += '<div id="exProjects">';
	dialog += "<p>Select the <strong>test data</strong>:</p>";
	dialog += '&nbsp;&nbsp;&nbsp;<input type="radio" name="data" value="groel" checked>';
	dialog += '&nbsp;&nbsp;' + getRefTestData("groel");
	dialog += '<br />';
	dialog += '&nbsp;&nbsp;&nbsp;<input type="radio" name="data" value="bpv">';
	dialog += '&nbsp;&nbsp;' + getRefTestData("bpv");
	dialog += '<br />';
	dialog += '&nbsp;&nbsp;&nbsp;<input type="radio" name="data" value="ribosome">';
	dialog += '&nbsp;&nbsp;' + getRefTestData("ribosome");
	dialog += '<br />';
	dialog += "</div>";
	dialog += "<br />";
	
	accessPopup(title, dialog, 'getProjExample', 'Go to project', 'Cancel');
		
//	dialog += '<div id="exProjects">';
//	dialog += "<p>Click over the <strong>example project</strong> to visit it:</p><br />";
//	dialog += '&nbsp;&nbsp;&nbsp;<a href="javascript:getProjExample(1);" style="color:firebrick;">·&nbsp;Groel data project</a>';
//	dialog += '&nbsp;&nbsp;' + getRefTestData("groel");
//	dialog += '<br />';
//	dialog += '&nbsp;&nbsp;&nbsp;<a href="javascript:getProjExample(2);" style="color:firebrick;">·&nbsp;BPV project</a>';
//	dialog += '&nbsp;&nbsp;' + getRefTestData("bpv");
//	dialog += '<br />';
//	dialog += '&nbsp;&nbsp;&nbsp;<a href="javascript:getProjExample(3);" style="color:firebrick;">·&nbsp;Ribosome project</a>';
//	dialog += '&nbsp;&nbsp;' + getRefTestData("ribosome");
//	dialog += '<br />';
//	dialog += "</div>";
	
//	new Messi(dialog, {
//		title : title,
//		modal : true,
//		buttons : [ {
//			id : 0,
//			label : 'Close',
//			val : 'X',
//			btnClass : 'fa-times'
//		} ]
//	});
}

function getProjExample(elm){
	var x = $("div#exProjects input[type='radio']:checked").val();
	switch(x){
		case "groel":
			var url = "/project_content/?projectName=GroelTestData&mode=service";
			break;
		case "bpv":
			var url ="/project_content/?projectName=BpvTestData&mode=service";
			break;
		case "ribosome":
			var url ="/project_content/?projectName=RiboTestData&mode=service";
			break;
	}
	goWithSubDomainURL(url);
}

function downloadTestdata(elm){
	var selected = $("#testData input[type='radio']:checked").val();
	var URL = getSubDomainURL() + "/get_testdata/?testData="+ selected;
	
	$.ajax({
		type : "GET",
		url : URL,
		dataType: "text",
		async: false,
		success : function(text) {
			var str = text.split("/")
			var fn = str[str.length-1]
			var URL = getSubDomainURL() + '/get_file/?path='+text+'&filename='+fn
			window.location.href = URL;
		}
	});
}

function createServProject(elm) {
	projName = randomString(32, '#aA')
	var selected = $("#testData input[type='radio']:checked").val();

	var URL = getSubDomainURL() + "/create_service_project/?projectName=" + projName
	if(selected != undefined){
		URL += "&testData="+selected;
	}
	
	$.ajax({
		type : "GET",
		url : URL,
		async: false,
		success : function() {
			var title = "ACCESS CODE"
			var msg = "<p>Your <strong>access code</strong> to the <strong>Project</strong> is:</p>" +
				"<br /><p><h2>" + projName + "</h2></p><br />" +
                "<p>Please <strong>save this code securely</strong> " +
				"to access the project in the future.</p>";
			var msg = msg + "<input type='hidden' class='content' value='" + projName + "' />";
			var funcName = "goToProject"

			accessPopup(title, msg, funcName, 'Go to the project', 'Exit');
		}
	});
}

function goToProjectForm() {
	var title = 'Confirm access code'
	var dialog = "<p>Please write the <strong>identification code</strong> " +
		"to access to your <strong>Project</strong>.</p>" +
		"<p>This code was given to you after the creation of your project.</p>" +
		"<p>If you forgot it, please contact the <a style='color:firebrick;' href='#'>Scipion group</a>.</p><br />";
	var msg = dialog + "<p><input type='text' id='code' class='content' style='width:100%;text-align:center;'/></p>";
	var funcName = 'goToProject';
	accessPopup(title, msg, funcName, 'Confirm', 'Cancel');
}

function goToProject(elm) {
	var code = elm.val();
	
	// remove the blank spaces
	code = code.replace(/\s+/g, '');
	
	var URL = getSubDomainURL() + "/check_project_id/?code=" + code;
	$.ajax({
		type : "GET",
		url : URL,
		success : function(result) {
			if (result == 1) {
				var URL2 = getSubDomainURL() + "/project_content/?projectName="+code+"&mode=service";
				window.location.href = URL2;
			} else {
				var title = "Bad Access";
				var msg = "<p>Wrong <strong>access code</strong>!</p>" +
					"<p>If you forgot the code, please contact the <a style='color:firebrick;' href='#'>Scipion group</a>.</p>";
				errorPopup(title, msg);
			}
		}
	});
}

function getRefTestData(id){
	var ref = ""
	switch(id){
		case "bpv":
			ref = "<strong>Bovine Papillomavirus</strong> (8 averages, 100x100 pixels, <a href='http://dx.doi.org/10.1073/pnas.0914604107' style='color:firebrick;' target='_blank'>from Wolf,M. et al. (2010)</a>)"
			break;
		case "groel":
			ref = "<strong>GroEL</strong> (44 averages, 64x64 pixels, <a href='http://dx.doi.org/10.1016/j.str.2004.05.006' style='color:firebrick;' target='_blank'>from Ludtke, S.J. et al. (2004)</a>)"
			break;
		case "ribosome":
			ref = "<strong>Eukaryotic Ribosome</strong> (32 averages, 64x64 pixels, <a href='ftp://ftp.ebi.ac.uk/pub/databases/emtest/SPIDER_FRANK_data/' style='color:firebrick;' target='_blank'>from J.Frank</a>)"
			break;
	}
	return ref;
}
