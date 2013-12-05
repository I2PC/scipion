<html>
<head>
  <title>Scipion Interactive Protocol Menu</title>
</head>
<body>
	<p>Opening interactive protocol:</p>
	<form name="interactiveProtocolMenuForm" action="/runAction" method="post">
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
				<td> Number of Trials</td>
				<td> <input type="text" name="numberTrials" value="{{numberTrials}}"></td>
			</tr>
			<!-- Hidden Fields -->
			<input type="hidden" name="pathAskpass" value="{{pathAskpass}}">
			<input type="hidden" name="action" value="runInteractiveProtocol">
			<tr>
				<td><input type="submit" value="Run Protocol"></td>
				<td><input type="button" value="Cancel"></td>
		</table>
	</form>

</body>
</html>


