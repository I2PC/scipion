<html>
<head>
  <title>Scipion Transfer Files Menu</title>
</head>
<body>
	<p>Opening Scipion Transfer Files Menu:</p>
	<form name="transferFilesMenuForm" action="/runTransferFiles" method="post">
		<table border="1">
			<tr>
				<td> Commands: </td>
				<td> <input type="text" name="commands" value="{{commands}}"></td>
			</tr>
			<tr>
				<td> Mode: </td>
				<td> <input type="text" name="commands" value="{{mode}}"></td>
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
			<tr>
				<td> Source Files</td>
				<td> <input type="text" name="filesSource" value="{{filesSource}}"></td>
			</tr>
			<tr>
				<td> Source Path</td>
				<td> <input type="text" name="sourcePath" value="{{sourcePath}}"></td>
			</tr>
			<tr>
				<td> Target Path</td>
				<td> <input type="text" name="targetPath" value="{{targetPath}}"></td>
			</tr>
			<tr>
				<td> Can Rsync</td>
				<td> <input type="text" name="canRsync" value="{{canRsync}}"></td>
			</tr>
			<tr>
				<td> Available Tunnel</td>
				<td> <input type="text" name="availableTunnel" value="{{availableTunnel}}"></td>
			</tr>
			<tr>
				<td> Port Tunnel</td>
				<td> <input type="text" name="portTunnel" value="{{portTunnel}}"></td>
			</tr>
			<!-- Hidden Fields -->
			<input type="hidden" name="pathAskpass" value="{{pathAskpass}}">
			<tr>
				<td><input type="submit" value="Run Protocol"></td>
				<td><input type="button" value="Cancel"></td>
		</table>
	</form>

</body>
</html>


