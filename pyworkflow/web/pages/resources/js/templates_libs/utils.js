/**
 * Generic lib with commons methods
 * 
 * popup(URL); 
 * customPopup(URL, widthValue, heightValue);
 * customPopupHTML(html);  
 * closePopup();
 * 
 */

function popup(URL) {
	var popup_width = 500;
	var popup_height = 470;
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

function customPopupHTML(html) {
	day = new Date();
	id = day.getTime();
	var popup = window.open('', id, 'height=470,width=825');
	popup.document.write(html);
}

function closePopup() {
	// opener.location.reload(true);
	// self.close();
	window.opener.location.reload(true);
	window.close();
}
