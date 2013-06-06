function createProjectForm() {

	var html = "Project Name: <input type='text' id='newProjName' class='content'/>";

	new Messi(html, {
		title : 'Enter the project name',
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Ok',
			val : 'Y',
			btnClass : 'btn-select',
			btnFunc : 'createProject'
		}, {
			id : 1,
			label : 'Cancel',
			val : 'C',
			btnClass : 'btn-cancel'
		} ]
	});
}

function createProject(elm) {
	var projName = elm.attr('value');
	$.ajax({
		type : "GET",
		url : "/create_project/?projectName=" + projName,
		success : function() {
			window.location.href = "/projects/";
		}
	});
}

function deleteProjectForm(projName) {

	var msg = "<div class='content' value='" + projName
			+ "'>Are you sure to <strong>DELETE</strong> project <strong>"
			+ projName + "</strong> and all its <strong>DATA</strong></div>";

	new Messi(msg, {
		title : 'Confirm project deletion',
		modal : true,
		buttons : [ {
			id : 0,
			label : 'Yes',
			val : 'Y',
			btnClass : 'btn-select',
			btnFunc : 'deleteProject'
		}, {
			id : 1,
			label : 'No',
			val : 'C',
			btnClass : 'btn-cancel'
		} ]
	});
}

function deleteProject(elm) {
	var projName = elm.attr('value');
	$.ajax({
		type : "GET",
		url : "/delete_project/?projectName=" + projName,
		success : function() {
			window.location.href = "/projects/";
		}
	});
}
