<head>
  
</head>

<?php
/********************************************
* printTable outputs any SQL table in an    *
* HTML tabular format.                      *
* ARGUMENTS:                                *
* $result :	the return value of a call  *
*		to mysql_query()            *
* RETURNS:                                  *
* true if there is data in the result set   *
* false if there is no data in the result   *
* set					    *
********************************************/
function printTable($result)
{
	//get the first row in the table:
	if(!$row = pg_fetch_assoc($result))
	{
		return false;
	}
	
	//set up the table:
	print("<table><tr>");
	
	//print the column names as table headers:
	foreach($row as $key=>$value)
	{
		print("<th>$key</th>");
	}
	print("</tr>");
	
	//loop through the table, printing
	//the field values in table cells:
	do
	{
		print("<tr>");
		foreach($row as $key=>$value)
		{
			print("<td>$value</td>");
		}
		print("</tr>");
	}
	while ($row = pg_fetch_assoc($result));
	
	//close out the table:
	print("</tr></table>");
	return true;
}
?>

<body>

<?php

require("../connectionStrings.php");

// create connection
$connection  = pg_connect($pgConnectStr1);
if (!$connection) {
echo "Couldn't make a connection!";
exit;
}

// test connection
if (!$connection) {
echo "Couldn't make a connection!";
exit;
}

// create SQL statement
$sql = "select * from xmipp_users order by my_date";

// execute SQL query and get result
//$sql_result = pg_exec($connection,$sql);
if(!$result = pg_query($sql))
{
   print ("Query could not be executed.\n");
}

$r = printTable($result); 
pg_close($connection);

?> 
</body>
