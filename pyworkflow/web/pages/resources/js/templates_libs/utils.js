 /*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
 * 			   Adrian Quintana (aquintana@cnb.csic.es)
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
 * Generic library with commons methods.
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
 * 
 * function detectWebBrowser()
 * 	->	Function to detect the web browser used in the navigation.
 * 
 * function popup(URL)
 * 	-> Launch a basic popup (600x500) opening the URL passed by argument.
 * 
 * function customPopup(URL, widthValue, heightValue)
 * 	->	Launch a popup opening the URL passed by argument. 
 * 		The size of the popup is customized with the width and height chosen.
 * 
 * function customPopupHTML(html, widthValue, heightValue)
 * 	->	Launch a popup with the HTML code passed by argument.
 *  	The size of the popup is customized with the width and height chosen.
 *  
 * function popUpJSON(json)
 * 	->	This method recive a JSON, and depending of the key content, open a 
 * 		diferent popups with diferent settings.
 * 		This function in the analyze results of the protocols runs.
 * 
 * function showPlot(url)
 * 	->	Function to show a xplotter(PNG) in a adjusted popup.
 * 
 * function getUrlParameters(parameter, staticURL, decode)
 * 	->	Auxiliar function to obtain individual parameters from a URL.
 * 
 * function infoPopup(title, msgText, autoclose, closeFunc) 
 * 	->	Creates a messi popup with a title and message passed by arguments.
 * 		If autoclose then the opener window will be closed when confirm button is pressed and closeFunc if provided will be executed.
 * 		It is used to show the help content in the protocol form and info message anywhere.
 * 
 * function warningPopup(title, msgText, funcName) 
 * 	->	Creates a messi popup with a title and message passed by arguments.
 * 		funcName is the function that will be executed if 'Yes' option is selected
 * 		It is used to show any warning message
 * 
 * function errorPopup(title, msgText) 
 * 	->	Creates a messi popup with a title and message passed by arguments.
 * 		It is used to show any error message
 * 
 * function isNaturalNumber(n)
 *  ->  Check if n is natural and returns true or false
 *
 * function editObjParam(id, title_label, value_label, title_comment, 
 * 			value_comment, msg_comment, typeObj)
 *  ->	Launch a messi popup with an input and textarea to edit the label and comment
 * 		for an object.
 * 
 * function updateLabelComment()
 * 	->	Method to store the label and comment for an object.
 * 
 *  * function launchViewer(id)
 * 	->	Launch the viewers to analyze the results for an object (protocol or object).
 * 
 ******************************************************************************/

/** METHODS *******************************************************************/

function startsWith(str, pattern){
	return str.lastIndexOf(pattern, 0) === 0
}


function detectWebBrowser(){
	/*
	 * Function to detect the web browser used in the navigation.
	 */
	
	var res = ""

	if(!!window.opera || navigator.userAgent.indexOf(' OPR/') >= 0){
	    // Opera 8.0+ (UA detection to detect Blink/v8-powered Opera)
		res = "opera";
	}
	else if(typeof InstallTrigger !== 'undefined'){
		// Firefox 1.0+
		res = "firefox";
	}
	else if(Object.prototype.toString.call(window.HTMLElement).indexOf('Constructor') > 0){
		// At least Safari 3+: "[object HTMLElementConstructor]"
		res = "safari"; 
	}
	else if(!!window.chrome && !isOpera){
		// Chrome 1+
		res = "chrome";
	}
	else if(/*@cc_on!@*/false || !!document.documentMode){
		// At least IE6
		res = "ie";
	}
	
	return res;
}


function popup(URL) {
	/*
	 * Launch a basic popup (600x500) opening the URL passed by argument.
	 */
	var popup_width = 600;
	var popup_height = 500;
	day = new Date();
	id = day.getTime();
	eval("page"
			+ id
			+ " = window.open(URL, '"
			+ id
			+ "', 'toolbar=0,scrollbars=1,location=0,statusbar=0,menubar=0,resizable=0,width='+popup_width+',height='+popup_height+'');");
}

function customPopup(URL, widthValue, heightValue) {
	/*
	 * Launch a popup opening the URL passed by argument. 
	 * The size of the popup is customized with the width and height chosen.
	 */
	day = new Date();
	id = day.getTime();
	eval("page"
			+ id
			+ " = window.open(URL, '"
			+ id
			+ "', 'toolbar=0,scrollbars=1,location=0,statusbar=0,menubar=0,resizable=0,width='+widthValue+',height='+heightValue+'');");
}

function customPopupFileHTML(html, title, widthValue, heightValue) {
	/*
	 * Launch a popup with the HTML code passed by argument.
	 * The size of the popup is customized with the width and height chosen.
	 */
	day = new Date();
	id = day.getTime();
	var popup = window.open('', id, 'height='+heightValue+',width='+widthValue);
	
	style = "background-color:black;color:white;font-family:Monospace;padding:1em;"
    title = "<title>"+ title + "</title>"
    body = "<div style="+ style +">" + title + html + "</div>"

	popup.document.write(body);
	popup.document.close();
}

function customPopupHTML(html, widthValue, heightValue) {
	/*
	 * Launch a popup with the HTML code passed by argument.
	 * The size of the popup is customized with the width and height chosen.
	 */
	day = new Date();
	id = day.getTime();
	var popup = window.open('', id, 'height='+heightValue+',width='+widthValue);
	popup.document.write(html);
	popup.document.close();
}

function customPopUpFile(url){
	$.ajax({
		type : "GET",
		url : url,
		dataType : "html",
		success : function(html) {
			customPopupHTML(html,600,500);
		}
	});	
}

function popUpJSON(json){
	/*
	 * This method receive a JSON, and depending of the key content, open a 
	 * diferent popups with diferent settings.
	 * This function in the analyze results of the protocols runs.
	 */

	// Open pop-ups depending of JSON parameters
	$.each(json, function(i, item) {
		array = item.split('::')
		key = array[0]
		value = array[1]
		//alert("item=" + item + " key=" + key + " value="+value)

		if (key=="url_form"){
			popup(value);
		} else if(key=="urls"){
			for(var x=0;x<value.length;x++){
				customPopup(value[x],1000,900);
			}
		} else if(key=="html"){
			customPopupHTML(value,600,500);
		} else if(key=="file"){
			customPopUpFile(value);
		} else if(key=="error"){
			errorPopup("Error",value);
		} else if(key=="showj"){
			customPopup(value,1024,800);
		} else {
			popup(value);
		}
	});
}

function getUrlParameters(parameter, staticURL, decode){
	/*
	 * Auxiliar function to obtain individual parameters from a URL.
	 */
   var currLocation = (staticURL.length)? staticURL : window.location.search,
       parArr = currLocation.split("?")[1].split("&"),
       returnBool = true;
   
   for(var i = 0; i < parArr.length; i++){
        parr = parArr[i].split("=");
        if(parr[0] == parameter){
            return (decode) ? decodeURIComponent(parr[1]) : parr[1];
            returnBool = true;
        }else{
            returnBool = false;            
        }
   }
   if(!returnBool) return false;  
}

function infoPopup(title, msgText, autoclose, closeFunc) {
	/*
	 * Creates a messi popup with a title and message passed by arguments.
	 * If autoclose then the opener window will be closed when confirm button is pressed and closeFunc if provided will be executed.
	 * It is used to show the help content in the protocol form and info message anywhere.
	 */
	
	//HTML to be used in the information popup
	msg="<table><tr><td><i class=\"fa fa-info-circle fa-4x\" style=\"color:#6fabb5;\"></i>"
		+ "</td><td>"+ msgText +"</td></tr></table>";
	
	
	if(autoclose){
		new Messi(msg, {
			title : title,
			modal : true,
			buttons : [ {
				id : 0,
				label : 'Close',
				val : 'X',
				btnClass : 'fa-times'
			} ],
			callback: function(){
				window.close();
				if (closeFunc){
					eval(closeFunc)
				}
			}
		});
	}
	else{
		new Messi(msg, {
			title : title,
			modal : true,
			buttons : [ {
				id : 0,
				label : 'Close',
				val : 'X',
				btnClass : 'fa-times'
			} ]
		});
	}	
}

function warningPopup(title, msgText, funcName){
	/*
	 * Creates a messi popup with a title and message passed by arguments.
	 * funcName is the function that will be executed if 'Yes' option is selected
	 * It is used to show any warning message
	 */
	
	// HTML to be used in the warning popup
	msg = "<table><tr><td><i class=\"fa fa-warning fa-4x\" style=\"color:#fad003;\"></i>"
		+ "</td><td>"+ msgText +"</td></tr></table>"

	new Messi(msg, {
		title : title,
		// modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : funcName
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'fa-ban'
		} ]
	});
}

function errorPopup(title, msgText){
	/*
	 * Creates a messi popup with a title and message passed by arguments.
	 * It is used to show any error message
	 */
	
	//HTML to be used in the error popup
	msg = "<table><tr><td><i class=\"fa fa-times-circle fa-4x\" style=\"color:firebrick;\"></i>"
	+ "</td><td>"+ msgText +"</td></tr></table>"
	
	new Messi(msg, {
		title : title,
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Close',
			val : 'X',
			btnClass : 'fa-times'
		} ]
	});
}

function listToString(list){
	var res = "";
	
	for (var x=0;x<list.length;x++){
		var elm = list[x];
		if(res.length == 0){
			res += elm;
		} else{
			res += "," + elm;
		}
	}
	
	if(res.length == 0){
		res="None";
	}
	
	return res;
}

function isNaturalNumber(n) {
    n = n.toString(); // force the value incase it is not
    var n1 = Math.abs(n),
        n2 = parseInt(n, 10);
    
    return !isNaN(n1) && n2 === n1 && n1.toString() === n && n2>0;
}


function processObjParam(id, title_label, title_comment,
                      msg_comment, typeObj) {
	
	// New values
	var value_label = $("input#runName").val();
	var value_comment= $("input#comment").val();
		
	editObjParam(id, title_label, value_label, 
	    title_comment, value_comment,
	    msg_comment, typeObj);

}

function editObjParam(id, title_label, value_label, 
                      title_comment, value_comment,
                      msg_comment, typeObj) {
	/*
	 * Launch a messi popup with an input and textarea to edit the label and comment
	 * for an object.
	 */
	
	// Create the form table to edit label and comment
	var html = "<table id='params' data-type='"+ typeObj +"'value='"+ id +"'>" + 
		 	"<tr>" +
		 	"<td>" +"<h3>"+ title_label +"</h3>" +"</td>" +
		 	"<td>" + "<input id='label_new' value='"+ value_label+"'>" + "</td>" +
			"</tr>"+
			"<tr>" +
		 	"<td>" +"<h3>" + title_comment +"</h3>" +"</td>" +
		 	"<td>" + "</h3><textarea id='comment_new'>"+ value_comment + "</textarea>" + "</td>" +
			"</tr>"+
			"</table>"
	
	// Some replaces to note up: &#013;&#010;
	
	// Launch the form with messi.js to start the edition
	new Messi(html, {
		title : 'Object Editor',
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Select',
			val : 'Y',
			btnClass : 'fa-check',
			btnFunc : "updateLabelComment"
			}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'fa-ban'
		}]
	});
}


function updateLabelComment(){
	/*
	 * Method to store the label and comment for an object.
	 */
	var elm_table = $("table#params");
	var id = elm_table.attr('value');
	var typeObj = elm_table.attr('data-type');
	
	// New values
	var value_label = $("input#label_new").val();
	var value_comment= $("textarea#comment_new").val();
		
	if (id == 'new'){
		// Return the new values to the source form
		returnLabelComment(value_label, value_comment);
		
	} else {
		var url_param = "/set_attributes/?" +
			"id=" + id + 
			"&label=" + value_label + 
			"&comment=" + value_comment 
			
		$.ajax({
			type : "GET",
			url : encodeURI(url_param),
			async: false,
			success : function() {
				// Return the new values to the source form
				returnLabelComment(value_label, value_comment);
			}
		});
		
	}
}

function returnLabelComment(label, comment){
	// Return the new values to the source form
	$("input#runName").val(label);
	$("input#comment").val(comment);
}


function replaceAll(find, replace, str) {
	return str.replace(new RegExp(find, 'g'), replace);
}

function launchViewer(id){
	/*
	 * Launch the viewers to analyze the results of the protocol run
	 */
	$.ajax({
		type : "GET",
		// Execute the viewer 
		url : "/launch_viewer/?objectId=" + id,
		dataType : "json",
		success : function(json) {
			popUpJSON(json);
		}
	});	
}
