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
 * Generic lib with commons methods.
 * 
 * ATTRIBUTES LIST:
 * 
 * METHODS LIST:
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
 * function closePopup()
 * 	->	Function to close the popup what call this method.
 * 
 * function showErrorValidation(json)
 * 	->	Function to normalize the errors launched in the protocol form when a 
 * 		protocol cannot be launched.
 * 
 * function launchMessiSimple(title, msg, autoclose)
 * 	->	Function to launch a basic simple window with a title and message passed
 * 		by argument. The autoclose can be chose with a boolean.
 * 
 * function messiError(msg)
 * 	->	This function return the HTML to be used in an error popup with messi.js
 * 
 * function messiWarning(msg)
 *  ->	This function return the HTML to be used in an warning popup with messi.js
 * 
 * function messiInfo(msg)
 *  ->	This function return the HTML to be used in an information popup with messi.js
 * 
 ******************************************************************************/

/** METHODS ******************************************************************/

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

function popUpJSON(json){
	/*
	 * This method recive a JSON, and depending of the key content, open a 
	 * diferent popups with diferent settings.
	 * This function in the analyze results of the protocols runs.
	 */
	// Open pop-ups depending of JSON parameters
	$.each(json, function(key, value) {
		if(key=="url_form"){
			customPopup(value,500,350);
		} else if(key=="showj"){
			customPopup(value,1024,600);
		} else if(key=="showjs" || key=="urls"){
			for(var x=0;x<value.length;x++){
				customPopup(value[x],1024,600);
			}
		} else if(key=="url"){
			customPopup(value,1024,600);
		} else if(key=="html"){
			customPopupHTML(value,600,500);
		} else if(key=="plots"){
			for(var x=0;x<value.length;x++){
				showPlot(value[x]);
			}
		} else if(key=="plotsComposite" || key=="plot"){
			showPlot(value);
		} else if(key=="error"){(
			launchMessiSimple("Error",messiError(value)));
		} else {
			customPopup(value,800,600);
		}
	});
}

function showPlot(url){
	/*
	 * Function to show a xplotter(PNG) in a adjusted popup.
	 */
	width = getUrlParameters("width", url, true)
	height = getUrlParameters("height", url, true)
//	alert(width +"x" + height)
	customPopup(url,width,height);
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


function closePopup() {
	/*
	 * Function to close the popup what call this method.
	 */
	// opener.location.reload(true);
	// self.close();
	//	window.opener.location.reload(true);
	window.close();
}

function showErrorValidation(json) {
	/*
	 * Function to normalize the errors launched in the protocol form when a 
	 * protocol cannot be launched.
	 */
	var msg = JSON.stringify(json);
	msg = msg.replace("<", "");
	msg = msg.replace(">", "");
	msg = msg.replace("[", "");
	msg = msg.replace("]", "");
	
	launchMessiSimple('Errors found',messiError(msg));
}

function launchMessiSimple(title, msg, autoclose){
	/*
	 * Function to launch a basic simple window with a title and message passed
	 * by argument. The autoclose can be chose with a boolean.
	 */
	if(autoclose){
		new Messi(msg, {
			title : title,
			modal : true,
			buttons : [ {
				id : 0,
				label : 'Ok',
				val : 'Y',
				btnClass : 'btn-select'
			}],
			callback: function(){
				closePopup();
			}
		});
	} else {
		new Messi(msg, {
			title : title,
			modal : true,
			buttons : [ {
				id : 0,
				label : 'Ok',
				val : 'Y',
				btnClass : 'btn-select'
			}]
		});
	}
}
	
function messiError(msg){
	/*
	 * This function return the HTML to be used in an error popup with messi.js
	 */
	var res = "<table><tr><td><i class=\"fa fa-times-circle fa-4x\" style=\"color:firebrick;\"></i>"
	+ "</td><td>"+ msg +"</td></tr></table>";

	return res;
}

function messiWarning(msg){
	/*
	 * This function return the HTML to be used in an warning popup with messi.js
	 */
	var res = "<table><tr><td><i class=\"fa fa-warning fa-4x\" style=\"color:#fad003;\"></i>"
	+ "</td><td>"+ msg +"</td></tr></table>";

	return res;
}

function messiInfo(msg){
	/*
	 * This function return the HTML to be used in an information popup with messi.js
	 */
	var res = "<table><tr><td><i class=\"fa fa-info-circle fa-4x\" style=\"color:#6fabb5;\"></i>"
	+ "</td><td>"+ msg +"</td></tr></table>";

	return res;
}
