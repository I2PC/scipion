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
 * function showErrorValidation(json)
 * 	->	Function to normalize the errors launched in the protocol form when a 
 * 		protocol cannot be launched.
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
		} else if(key=="error"){
				errorPopup("Error",value);
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
	
	errorPopup('Errors found',msg);
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
			title : 'Help' + ' ' + title,
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
			title : 'Help' + ' ' + title,
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
	
	//HTML to be used in the warning popup
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
