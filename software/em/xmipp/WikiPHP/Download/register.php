<?php

$errmsg       = ''; // error message
$sname        = ''; // sender's name
$email        = ''; // sender's email address
$organization = ''; //
$mailoption   = ''; // mailoption
$workshop   = ''; // workshop to be registered in
//$workshop   	 = $_POST['workshop'];

//print $_POST['workshop'];
parse_str($_SERVER['QUERY_STRING']);
//print $_SERVER['REQUEST_URI'];
//print "$workshop<br/>";


//if(isset($_POST['sname']) || isset($_POST['email']) || isset($_POST['organization'])){ 
//if(($_POST['sname']!='') || ($_POST['email']!='') || ($_POST['organization']!='')){ 

if(isset($_POST['send'])){
  $sname         = $_POST['sname'];
  $email         = $_POST['email'];
  $organization  = $_POST['organization'];
  $workshop   	 = $_POST['workshop'];
  $direction     = $_SERVER['REMOTE_ADDR'];
  $mailoption    = $_POST['mailoption'];

  print "$sname<br/>";
  print "$email<br/>";
  print "$organization<br/>";
  print "$workshop<br/>";

  //if((trim($sname) != '') || (trim($organization) != '') || (trim($email) != '')) {
  if(trim($sname) == ''){
    $errmsg = 'Please enter your name';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/RegisterForm?workshop=' . "$workshop" . '&errormsg=' . "$errmsg" ) ;
  } else if(trim($organization) == '') {
    $errmsg = 'Please enter your organization';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/RegisterForm?workshop=' . "$workshop" . '&errormsg=' . "$errmsg" ) ;
  } else if(trim($email) == ''){
    $errmsg = 'Please enter your email address';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/RegisterForm?workshop=' . "$workshop" . '&errormsg=' . "$errmsg" ) ;
  } else if(!isEmail($email)){
    $errmsg = 'Your email address is not valid';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/RegisterForm?workshop=' . "$workshop" . '&errormsg=' . "$errmsg" ) ;
  } else if(trim($mailoption) == ''){
    $errmsg = 'Please select an email option';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/RegisterForm?workshop=' . "$workshop" . '&errormsg=' . "$errmsg" ) ;
  }

  if($errmsg == ''){
    if(get_magic_quotes_gpc()){
      $sname        = stripslashes($sname);
      $email        = stripslashes($email);
      $organization = stripslashes($organization);
      $mailoption   = stripslashes($mailoption);
      $direction    = stripslashes($direction);
      $workshop	    = stripslashes($workshop);
    }	
    $sname        = pg_escape_string($sname);
    $email        = pg_escape_string($email);
    $organization = pg_escape_string($organization);
    $mailoption   = pg_escape_string($mailoption);
    $direction    = pg_escape_string($direction);
    $workshop	  = pg_escape_string($workshop);

    // extra validations

    $date = date('m/d/Y h:i:s a', time());
    //Data insertion into the userTable table (xmipp database)
    $sql_query = "INSERT INTO userTable (sname,email,organization,workshop, direction, timestamp, mailoption) 
    $conn = pg_connect("host=xxxx port=xxxx dbname=xxxx user=xxxxxxn password=xxxxxxx");

    VALUES('$sname', '$email','$organization','$workshop', '$direction', '$date', '$mailoption');";
    $res = pg_query($sql_query);

    if ( FALSE === $res ){
      print "Data entry error: " . pg_last_error() . "<br />";
      print "Contact xmipp@cnb.uam.es". "<br />";
      print $sql_query . "<br />";
    }

    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/registerConfirmated?workshop=' . "$workshop") ;
  } 
//}
}


//if(!isset($_POST['send']) || $errmsg != ''){
if((trim($sname) == '') || (trim($organization) == '') || (trim($email) == '') || $errmsg != ''){
?>

<form  action="register.php" method="post" name="msgform" id="msgform">
  <table width="500" border="0" align="center" cellpadding="2" cellspacing="1" class="twikiTable">
    <tr class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0"> 
      <td width="306" bgcolor="#ffffff" valign="top" class="twikiTableCol0 twikiFirstCol">Your Full Name</td>
      <td width="381" bgcolor="#ffffff" valign="top" class="twikiTableCol1 twikiLastCol"><input name="sname" type="text" class="box" id="sname" size="30" ></td>
    </tr>
    <tr class="twikiTableEven twikiTableRowdataBgSorted1 twikiTableRowdataBg1">
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol0 twikiFirstCol">Your Organization</td>
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol1 twikiLastCol"><input name="organization" type="text" class="box" id="organization" size="30"></td>
    </tr>
    <tr class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0"> 
      <td bgcolor="#ffffff" valign="top" class="twikiTableCol0 twikiFirstCol">Your Email</td>
      <td bgcolor="#ffffff" valign="top" class="twikiTableCol1 twikiLastCol"><input name="email" type="text" class="box" id="email" size="30"> <br>
    </tr>
    <tr class="twikiTableEven twikiTableRowdataBgSorted1 twikiTableRowdataBg1"> 
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol0 twikiFirstCol"></td>
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol1 twikiLastCol">
      <input type="radio" name="mailoption" value="nomail"> Don't subscribe me to news mail list<br>
      <input type="radio" name="mailoption" value="mail" checked> Subscribe me to news mail list<br>
      </td>
    </tr>
  <tr align="center" class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0">
    <td colspan="2"><input name="send" type="submit" class="twikiSubmit" id="send" value="Register" onclick="return checkForm();"></td>
  </tr>
  </table>
  <input type = "hidden" name="workshop" value="<?=$workshop ?>">
</form>
If you need any help you may contact us on <a href="mailto:xmipp@cnb.csic.es"> xmipp@cnb.csic.es </a>
<?php
}

function isEmail($email){
  return(preg_match("/^[-_.[:alnum:]]+@((([[:alnum:]]|[[:alnum:]][[:alnum:]-]*[[:alnum:]])\.)+(ad|ae|aero|af|ag|ai|al|am|an|ao|aq|ar|arpa|as|at|au|aw|az|ba|bb|bd|be|bf|bg|bh|bi|biz|bj|bm|bn|bo|br|bs|bt|bv|bw|by|bz|ca|cc|cd|cf|cg|ch|ci|ck|cl|cm|cn|co|com|coop|cr|cs|cu|cv|cx|cy|cz|de|dj|dk|dm|do|dz|ec|edu|ee|eg|eh|er|es|et|eu|fi|fj|fk|fm|fo|fr|ga|gb|gd|ge|gf|gh|gi|gl|gm|gn|gov|gp|gq|gr|gs|gt|gu|gw|gy|hk|hm|hn|hr|ht|hu|id|ie|il|in|info|int|io|iq|ir|is|it|jm|jo|jp|ke|kg|kh|ki|km|kn|kp|kr|kw|ky|kz|la|lb|lc|li|lk|lr|ls|lt|lu|lv|ly|ma|mc|md|mg|mh|mil|mk|ml|mm|mn|mo|mp|mq|mr|ms|mt|mu|museum|mv|mw|mx|my|mz|na|name|nc|ne|net|nf|ng|ni|nl|no|np|nr|nt|nu|nz|om|org|pa|pe|pf|pg|ph|pk|pl|pm|pn|pr|pro|ps|pt|pw|py|qa|re|ro|ru|rw|sa|sb|sc|sd|se|sg|sh|si|sj|sk|sl|sm|sn|so|sr|st|su|sv|sy|sz|tc|td|tf|tg|th|tj|tk|tm|tn|to|tp|tr|tt|tv|tw|tz|ua|ug|uk|um|us|uy|uz|va|vc|ve|vg|vi|vn|vu|wf|ws|ye|yt|yu|za|zm|zw)$|(([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5])\.){3}([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5]))$/i",$email));
}
?>

