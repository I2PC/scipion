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
function printTable($result,$Nrows)
{
	//get the first row in the table:
	if(!$row = pg_fetch_assoc($result))
	{
		return false;
	}
	
	//set up the table:
	print("<wholething><Nrows>$Nrows</Nrows>\n");
	
	//loop through the table, printing
	//the field values in table cells:
	do
	{
		print("<mrow>\n");
		foreach($row as $key=>$value)
		{
			print("<$key>$value</$key>\n");
		}
		print("</mrow>\n\n");
	}
	while ($row = pg_fetch_assoc($result));
        print("</wholething>");	
	return true;
	//close out the table:
}
?>
<?php
// create connection
require("../connectionStrings.php");
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
$sql = "select * from xmipp_users order by my_date limit 3";

// execute SQL query and get result
//$sql_result = pg_exec($connection,$sql);
if(!$result = pg_query($sql))
{
   print ("Query could not be executed.\n");
}
$Nrows = pg_num_rows($result);

$r = printTable($result,$Nrows); 
pg_close($connection);

?> 
