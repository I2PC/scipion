function launchToolbar(projName, id, elm) {
	var row = $("div#toolbar");

	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("tr#" + row.attr('value'));
		rowOld.attr('style', 'background-color: #fafafa;');
	}
	row.attr('value', id);
	$("tr#" + id).attr('style', 'background-color: #f2f2f2;');

	// Action Edit Button
	$("a#editTool").attr(
			'href',
			'javascript:popup("/form/?projectName=' + projName + '&protocolId='
					+ id + '")');
	// Action Copy Button
	$("a#copyTool").attr(
			'href',
			'javascript:popup("/form/?projectName=' + projName + '&protocolId='
					+ id + '")');
	// Action Copy Button
	$("a#deleteTool").attr(
			'href',
			'javascript:popup("/form/?projectName=' + projName + '&protocolId='
					+ id + '")');
	// Action Browse Button
	$("a#browseTool").attr(
			'href',
			'javascript:popup("/form/?projectName=' + projName + '&protocolId='
					+ id + '")');

	row.show(); // Show toolbar
}

function selTableMessi(elm) {

	var row = $("table.content");
	var id = elm.attr('id');

	if (row.attr('value') != undefined && row.attr('value') != id) {
		var rowOld = $("td#" + row.attr('value'));
		rowOld.attr('style', '');
	}
	row.attr('value', id);
	elm.attr('style', 'background-color: firebrick; color:white;');
}
