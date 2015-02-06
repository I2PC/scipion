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

function serviceProjForm(){
	var title = 'Project creation'
	var dialog = "<p>Your <strong>Project</strong> will be created.<br /><br />" +
        "This process generates a unique <strong>url access</strong>.<br /><br />" +
        "This url access should be used to have access to your data in future sessions.</p>" +
        "<p><br /></p>";
	
    dialog += "<p>Confirm to generate it.</p>";

	var funcName = 'createServProject';

	accessPopup(title, dialog, funcName, 'Confirm', 'Cancel');
}

function serviceTestDataForm(){
	var title = 'Test data'
	var dialog = ""
		
	dialog += "<p>You have two options to use <strong>Test data</strong> :<br />" +
        "1.- <strong>Create a project</strong> with <strong>Test data</strong> already imported inside.<br />" +	
		"2.- <strong>Download</strong> test files to your computer, to be manually imported into an already created project.</p>";
	dialog += "<br />";
	dialog += '<div id="testData">';
	dialog += "<p>Select <strong>Test data</strong>:</p>";
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
    dialog += "<p>After selection, choose your option.</p>";

    
    var btn1 = 'Create project'
	var ico1 = 'fa-check'
	var funcName1 = 'createServProject';	
		
	var btn2 = 'Download'
	var ico2 = 'fa-download';
	var funcName2 = 'downloadTestdata';
	
	accessPopup2opt(title, dialog, 
					 btn1, ico1, funcName1, 
					 btn2, ico2, funcName2, 
					 "Cancel")
}

function goExampleForm(){
	var title = 'Example projects'
	var dialog = ""

	dialog += '<div id="exProjects">';
	dialog += "<p>Select the <strong>Test data</strong>:</p>";
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
		
}

function getProjExample(elm){
	var x = $("div#exProjects input[type='radio']:checked").val();
	switch(x){
		case "groel":
			var url = "/service_content/?p=GroelTestData";
			break;
		case "bpv":
			var url ="/service_content/?p=BpvTestData";
			break;
		case "ribosome":
			var url ="/service_content/?p=RiboTestData";
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
			var funcHref = "javascript:goToProject(jQuery('.content'));"
				
			var msg = "<p>Your <strong>url to access </strong> to this <strong>Project</strong> is:</p>" +
			"<br /><p><h3>" +
			
			"<a style='color:firebrick;' href='http://scipion.cnb.csic.es/myfirstmap/service_content/?p="+ projName+ "'>" +
			
			"<a style='color:firebrick;' href='http://scipion.cnb.csic.es/myfirstmap/service_content/?p="+ projName+ "'>" +
			"http://scipion.cnb.csic.es/myfirstmap/service_content/?p="+ projName+ "</a>"+
			"</h3></p><br />" +
            "<p>Please <strong>save this url securely</strong> " +
			"in order to access to this project in future sessions.</p><br />";
			
			msg = msg + "<input type='hidden' class='content' value='" + projName + "' />";
			var funcName = "goToProject"

			accessPopup(title, msg, funcName, 'Go to the project', 'Exit');
		}
	});
}

//function goToProjectForm() {
//	var title = 'Confirm access code'
//	var dialog = "<p>Please write the <strong>identification code</strong> " +
//		"to access to your <strong>Project</strong>.</p>" +
//		"<p>This code was given to you after the creation of your project.</p>" +
//		"<p>If you forgot it, please contact us using this mail: <span style='color:firebrick;'>myfirstmap@cnb.csic.es</span></p><br />";
//	var msg = dialog + "<p><input type='text' id='code' class='content' style='width:100%;text-align:center;'/></p>";
//	var funcName = 'goToProject';
//	accessPopup(title, msg, funcName, 'Confirm', 'Cancel');
//}

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
				var URL2 = getSubDomainURL() + "/service_content/?p="+code;
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
			ref = "<strong>Eukaryotic Ribosome</strong> (32 averages, 64x64 pixels, <a href='ftp://ftp.ebi.ac.uk/pub/databases/emtest/SPIDER_FRANK_data/' style='color:firebrick;' target='_blank'>from J.Frank lab</a>)"
			break;
	}
	return ref;
}
