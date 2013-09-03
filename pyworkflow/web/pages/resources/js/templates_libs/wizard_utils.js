function selectList(elm) {
	var row = $("table#list");
	var oldValue = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != oldValue) {
		// unmark the last option
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', 'no-selected');
	}
	// mark the new option
	row.attr('value', oldValue);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');
	
	// load the image 
	var uri = elm.attr('value');
	if(uri == undefined){
		uri = elm.val();
	}
	
	
	
	
	
}
