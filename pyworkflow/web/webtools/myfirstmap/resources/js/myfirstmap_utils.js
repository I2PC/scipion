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

function serviceProjForm(){
	var title = 'Project creation'
	var dialog = "<p>Your <strong>project</strong> will be created with a unique <strong>url access</strong>.<br />" +
        "Use this URL to access your data in future sessions.</p>" +
        "<p><br /></p>";
	
    dialog += "<p>Confirm to generate it.</p>";

	var funcName = 'createServProject';

	accessPopup(title, dialog, funcName, 'Confirm', 'Cancel');
}

function serviceTestDataForm(){
	var title = 'Test data'
	var dialog = ""
		
	dialog += "<strong>Create a project</strong> with <strong>Test data</strong> already imported inside.<br /></p>";
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

    
    var btn1 = 'Create project'
	var ico1 = 'fa-check'
	var funcName1 = 'createServProject';	
		
		
	accessPopup(title, dialog, funcName1, btn1, "Cancel")
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
			var url = "/content/?p=GroelTestData";
			break;
		case "bpv":
			var url ="/content/?p=BpvTestData";
			break;
		case "ribosome":
			var url ="/content/?p=RiboTestData";
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
	var projName = "map"+randomString(16, '#aA')
	var selected = $("#testData input[type='radio']:checked").val();
	
	var projectUrl = "http://" + document.domain + getSubDomainURL() + "/m_content/?p="+ projName
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
			
			var msg = "<p>Your <strong>url to access </strong> this <strong>project</strong> is:</p>" +
			"<br /><p><h3>" + 
			"<a style='color:firebrick;' href='"+ projectUrl + "'>" +
			projectUrl + "</a>"+
			"</h3></p><br />" +
			"<p>The project will be <strong>DELETED TWO WEEKS</strong> after its creation.</p><br />"+
            "<p>Please <strong>SAVE or BOOKMARK this url securely</strong> " +
			"in order to access project in future sessions.</p>"+
			"<p>If you experience any problem contact us on this email: <span style='color:firebrick;'>scipion at cnb.csic.es</span></p>";
			
			msg = msg + "<input type='hidden' class='content' value='" + projName + "' />";
			var funcName = "goToProject"

			accessPopup(title, msg, funcName, 'Go to the project', 'Exit');
		}
	});
}

function goToProject(elm) {
	var code = elm.val();
	
	var URL2 = getSubDomainURL() + "/content/?p="+code;
	window.location.href = URL2;
	
	/*
	
	// remove the blank spaces
	code = code.replace(/\s+/g, '');
	
	var URL = getSubDomainURL() + "/check_project_id/?code=" + code+"&service=myfirstmap";
	$.ajax({
		type : "GET",
		url : URL,
		success : function(result) {
			if (result == 1) {
				var URL2 = getSubDomainURL() + "/content/?p="+code;
				window.location.href = URL2;
			} else {
				var title = "Bad Access";
				var msg = "<p>Wrong <strong>access code</strong>!</p>" +
					"<p>If you forgot the code, please contact the <a style='color:firebrick;' href='#'>Scipion group</a>.</p>";
				errorPopup(title, msg);
			}
		}
	});
	*/
}

function getRefTestData(id){
	var ref = ""
	switch(id){
		case "bpv":
			ref = "<strong>Bovine Papillomavirus</strong> (8 averages, 100x100 pixels, <a href='http://dx.doi.org/10.1073/pnas.0914604107' style='color:firebrick;' target='_blank'>from Wolf,M. et al. PNAS, 2010.</a>)"
			break;
		case "groel":
			ref = "<strong>GroEL</strong> (44 averages, 64x64 pixels, <a href='http://dx.doi.org/10.1016/j.str.2004.05.006' style='color:firebrick;' target='_blank'>from Ludtke, S.J. et al. Structure, 2004</a>)"
			break;
		case "ribosome":
			ref = "<strong>Eukaryotic Ribosome</strong> (32 averages, 64x64 pixels, <a href='ftp://ftp.ebi.ac.uk/pub/databases/emtest/SPIDER_FRANK_data/' style='color:firebrick;' target='_blank'>from J.Frank lab, data at EBI.</a>)"
			break;
	}
	return ref;
}
