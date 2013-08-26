function selectList(elm) {
	var row = $("div.list");
	var oldValue = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != oldValue) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
		rowOld.attr('class', 'no-selected');
	}
	row.attr('value', oldValue);
	elm.attr('style', 'background-color: LightSteelBlue;');
	elm.attr('class', 'selected');
	
	
}
