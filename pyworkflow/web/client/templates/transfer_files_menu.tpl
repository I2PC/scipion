<html>
<head>
  <title>Scipion Interactive Protocol Menu</title>
</head>
<body>
	<p>Opening interactive protocol:</p>
	<form name="interactiveProtocolMenuForm" action="/runInteractiveProtocol" method="post">
		<table border="1">
			<tr>
				<td> Commands: </td>
				<td> <input type="text" name="commands" value="{{commands}}"></td>
			</tr>
			<tr>
				<td> Machine: </td>
				<td> <input type="text" name="machine" value="{{machine}}"></td>
			</tr>
			<tr>
				<td> User: </td>
				<td> <input type="text" name="user" value="{{user}}"></td>
			</tr>
			<tr>
				<td> Password</td>
				<td> <input type="password" name="password" value="{{password}}"></td>
			</tr>
			<tr>
				<td><input type="submit" value="Run Protocol"></td>
				<td><input type="button" value="Cancel"></td>
		</table>
	</form>

</body>
</html>


