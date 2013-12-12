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
 * Generic lib with commons methods
 * 
 * popup(URL); 
 * customPopup(URL, widthValue, heightValue);
 * customPopupHTML(html);  
 * closePopup();
 * launchMessiSimple(title, msg);
 * messiError(msg);
 * messiWarning(msg);
 * messiInfo(msg);
 * 
 */

function popup(URL) {
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
	day = new Date();
	id = day.getTime();
	eval("page"
			+ id
			+ " = window.open(URL, '"
			+ id
			+ "', 'toolbar=0,scrollbars=1,location=0,statusbar=0,menubar=0,resizable=0,width='+widthValue+',height='+heightValue+'');");
}

function customPopupHTML(html, widthValue, heightValue) {
	day = new Date();
	id = day.getTime();
	var popup = window.open('', id, 'height='+heightValue+',width='+widthValue);
	popup.document.write(html);
	popup.document.close();
	
}

function popUpJSON(json){
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
		}else if(key=="url"){
			customPopup(value,1024,600);
		} else if(key=="html"){
			customPopupHTML(value,600,500);
		} else if(key=="plot"){
			customPopup(value,600,500);
		} else if(key=="plots"){
			for(var x=0;x<value.length;x++){
				customPopup(value[x],600,500)
			}
		} else if(key=="error"){
			launchMessiSimple("Error",messiError(value));
		} else {
			customPopup(value,800,600);
		}
	});
}


function closePopup() {
	// opener.location.reload(true);
	// self.close();
	//	window.opener.location.reload(true);
	window.close();
}

function showErrorValidation(json) {
	var msg = JSON.stringify(json);
	msg = msg.replace("<", "");
	msg = msg.replace(">", "");
	msg = msg.replace("[", "");
	msg = msg.replace("]", "");
	
	launchMessiSimple('Errors found',messiError(msg));
}

function launchMessiSimple(title, msg, autoclose){
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
				window.close();
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
	var res = "<table><tr><td><i class=\"fa fa-times-circle fa-4x\" style=\"color:firebrick;\"></i>"
	+ "</td><td>"+ msg +"</td></tr></table>";

	return res;
}

function messiWarning(msg){
	var res = "<table><tr><td><i class=\"fa fa-warning fa-4x\" style=\"color:#fad003;\"></i>"
	+ "</td><td>"+ msg +"</td></tr></table>";

	return res;
}

function messiInfo(msg){
	var res = "<table><tr><td><i class=\"fa fa-info-circle fa-4x\" style=\"color:#6fabb5;\"></i>"
	+ "</td><td>"+ msg +"</td></tr></table>";

	return res;
}
