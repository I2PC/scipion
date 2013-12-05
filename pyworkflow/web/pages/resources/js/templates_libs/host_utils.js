/*
 * Toolbar used in the view host template
 */
function launchHostsToolbar(projName, hostId, elm) {
	var row = $("#toolbar");

	if (row.attr('value') != undefined && row.attr('value') != hostId) {
		var rowOld = $("tr#" + row.attr('value'));
//		rowOld.attr('style', 'background-color: #fafafa;');
//		rowOld.attr('class', 'runtr');
		rowOld.removeClass('selected')
	}
	row.attr('value', hostId);
//	elm.attr('style', 'background-color: LightSteelBlue;');
//	elm.attr('class', 'selected');
	elm.addClass('selected')

	// Action Edit Button
	$("a#editTool").attr('href', 'javascript:editHost()');
	// Action Copy Button
	$("a#newTool").attr('href', 'javascript:newHost()');
	// Action Delete Button
	$("a#deleteTool").attr('href', 'javascript:deleteHost()');

	row.show(); // Show toolbar
}

	/* {% if message != None and message != '' %}
		window.onload = function() {
			alert("{{message}}")
		}
		{% endif %} */

var selectedHost = ""

function newHost() {
	document.location = "javascript:popup('/host_form/')"
}

function selectHost(hostId) {
	selectedHost = hostId
	document.getElementById("editTool").disabled = false;
	document.getElementById("deleteTool").disabled = false;
}

function editHost() {
	document.location = "javascript:popup('/host_form/?hostId=" + selectedHost
			+ "')";
}

function deleteHost() {
	if (selectedHost == "") {
		alert("No host has been selected to delete");
	} else {
		if (confirm("Are you sure you want to remove this host?")) {
			document.location.href = "/delete_host/?hostId=" + selectedHost;

		} else {
			alert("Remove operation cancelled");
		}
	}
}

/*
 * function getHost(hostLabel, projectName) { $.ajax({ type : "GET", url :
 * "/getHost/?hostLabel=" + hostLabel, dataType : "json", success :
 * function(json) { // specifying a dataType of json makes jQuery pre-eval the
 * response // for us ExecutionHostConfig host = json.host;
 * document.getElementById("scpnHosts").disabled = true;
 * document.getElementById("label").value = host.label;
 * document.getElementById("label").disabled = true;
 * document.getElementById("hostName").value = host.hostName;
 * document.getElementById("userName").value = host.userName;
 * document.getElementById("hostPath").value = host.hostPath; } }); }
 */
function clearForm() {
	document.getElementById("scpnHosts").disabled = false;
	document.getElementById("label").value = ""
	document.getElementById("label").disabled = false;
	document.getElementById("hostName").value = ""
	document.getElementById("userName").value = ""
	document.getElementById("hostPath").value = ""
}
function changeScpnHostSelection() {
	if (document.getElementById("scpnHosts").value != "") {
		alert("Not implemented yet")
	}
}