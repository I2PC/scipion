<html>
<head>
<TITLE>Xmipp Download Form</TITLE>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
<style type="text/css">
body {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 12px;
}

.box {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 12px;
	border: 1px solid #000000;

}

.bluebox {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 12px;
	font-weight: bolder;
	color: #FFFFFF;
	background-color: #006699;
	border: 1px solid #000000;
}

.maincell {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 12px;
	padding: 5px;
	border: 1px solid #006699;
}

.note {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 12px;
	padding: 5px;
	color: #CC0000;
}

.bigmaincell {
	font-family: Arial, Helvetica, sans-serif;
	font-size: 16px;
	padding: 5px;
	border: 1px solid #006699;
}

.errmsg {
	font-family: "Courier New", Courier, mono;
	font-size: 12px;
	font-weight: bolder;
	color: #CC0000;
}

</style>
<script language="JavaScript">
function checkForm()
{
	var cname, cemail, corganization, ccountry,cplatform, cversion, cmailoption;
	with(window.document.msgform)
	{
		cname      = sname;
		cemail     = email;
		cplatform  = platform;
		cversion   = version;
		ccountry   = country;
		corganization = organization;
		cmailoption = mailoption;
	}
	
	if(trim(cname.value) == '')
	{
		alert('Please enter your name');
		cname.focus();
		return false;
	}
	else if(trim(cemail.value) == '')
	{
		alert('Please enter your email');
		cemail.focus();
		return false;
	}
	else if(!isEmail(trim(cemail.value)))
	{
		alert('Email address is not valid');
		cemail.focus();
		return false;
	}
	else if(trim(corganization.value) == '')
	{
		alert('Please enter your organization');
		corganization.focus();
		return false;
	}
	else if(trim(ccountry.value) == '')
	{
		alert('Please select a country');
		ccountry.focus();
		return false;
	}
	else if(trim(cversion.value) == '')
	{
		alert('Please select xmipp version');
		cversion.focus();
		return false;
	}
	else if(trim(cplatform.value) == '')
	{
		alert('Please choose platform');
		cplatform.focus();
		return false;
	}
	else
	{
		cname.value          = trim(cname.value);
		cemail.value         = trim(cemail.value);
		corganization.value  = trim(corganization.value);
		cplatform.value      = trim(cplatform.value);
		cversion.value       = trim(cversion.value);
		ccountry.value       = trim(ccountry.value);
		cmailoption          = trim(cmailoption.value);
		return true;
	}
}

/*
Strip whitespace from the beginning and end of a string
Input : a string
*/
function trim(str)
{
	return str.replace(/^\s+|\s+$/g,'');
}

/*
Check if a string is in valid email format. 
Returns true if valid, false otherwise.
*/
function isEmail(str)
{
	var regex = /^[-_.a-z0-9]+@(([-_a-z0-9]+\.)+(ad|ae|aero|af|ag|ai|al|am|an|ao|aq|ar|arpa|as|at|au|aw|az|ba|bb|bd|be|bf|bg|bh|bi|biz|bj|bm|bn|bo|br|bs|bt|bv|bw|by|bz|ca|cat|cc|cd|cf|cg|ch|ci|ck|cl|cm|cn|co|com|coop|cr|cs|cu|cv|cx|cy|cz|de|dj|dk|dm|do|dz|ec|edu|ee|eg|eh|er|es|et|eu|fi|fj|fk|fm|fo|fr|ga|gb|gd|ge|gf|gh|gi|gl|gm|gn|gov|gp|gq|gr|gs|gt|gu|gw|gy|hk|hm|hn|hr|ht|hu|id|ie|il|in|info|int|io|iq|ir|is|it|jm|jo|jp|ke|kg|kh|ki|km|kn|kp|kr|kw|ky|kz|la|lb|lc|li|lk|lr|ls|lt|lu|lv|ly|ma|mc|md|mg|mh|mil|mk|ml|mm|mn|mo|mp|mq|mr|ms|mt|mu|museum|mv|mw|mx|my|mz|na|name|nc|ne|net|nf|ng|ni|nl|no|np|nr|nt|nu|nz|om|org|pa|pe|pf|pg|ph|pk|pl|pm|pn|pr|pro|ps|pt|pw|py|qa|re|ro|ru|rw|sa|sb|sc|sd|se|sg|sh|si|sj|sk|sl|sm|sn|so|sr|st|su|sv|sy|sz|tc|td|tf|tg|th|tj|tk|tm|tn|to|tp|tr|tt|tv|tw|tz|ua|ug|uk|um|us|uy|uz|va|vc|ve|vg|vi|vn|vu|wf|ws|ye|yt|yu|za|zm|zw)|(([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5])\.){3}([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5]))$/i;
	return regex.test(str);
}
</script>
</head>

<body>
<?php

$errmsg       = ''; // error message
$sname        = ''; // sender's name
$email        = ''; // sender's email addres
$organization = ''; //
$platform     = ''; // binaries compiled for... platform
$version      = ''; // xmipp version
$country      = ''; // country
$mailoption   = ''; // mailoption

if(isset($_POST['send']))
{
	$sname         = $_POST['sname'];
	$email         = $_POST['email'];
	$organization  = $_POST['organization'];
	$platform      = $_POST['platform'];
	$version       = $_POST['version'];
	$country       = $_POST['country'];
	$mailoption    = $_POST['mailoption'];
        $direction     = $_SERVER['REMOTE_ADDR'];
	
	require("connectionStrings.php");
        $conn = pg_connect($pgConnectStr2);
	
	if(trim($sname) == '')
	{
		$errmsg = 'Please enter your name';
	} 
	else if(trim($email) == '')
	{
		$errmsg = 'Please enter your email address';
	}
	else if(!isEmail($email))
	{
		$errmsg = 'Your email address is not valid';
	}
	else if(trim($platform) == '')
	{
		$errmsg = 'Please select platform';
	}
	else if(trim($organization) == '')
	{
		$errmsg = 'Please enter your organization';
	}
	else if(trim($country) == '')
	{
		$errmsg = 'Please select your country';
	}
	else if(trim($version) == '')
	{
		$errmsg = 'Please select Xmipp version';
	}
	else if(trim($mailoption) == '')
	{
		$errmsg = 'Please select an email option';
	}
	
	if($errmsg == '')
	{
		if(get_magic_quotes_gpc())
		{
			$sname        = stripslashes($sname);	  
			$email        = stripslashes($email);	  
			$organization = stripslashes($organization);
			$platform     = stripslashes($platform);    
			$version      = stripslashes($version);	  
			$country      = stripslashes($country);	  
			$mailoption   = stripslashes($mailoption);	  
		        $direction    = stripslashes($direction);	
		}	
	        $sname        = pg_escape_string($sname);      
		$email        = pg_escape_string($email);       
		$organization = pg_escape_string($organization);
		$platform     = pg_escape_string($platform);    
		$version      = pg_escape_string($version);     
		$country      = pg_escape_string($country);        
		$mailoption   = pg_escape_string($mailoption);        
		$direction    = pg_escape_string($direction);        
        	$sql_query = "INSERT INTO xmipp_users (sname, email,organization,
	                        	   platform,version,country,mailoption,direction) 
		VALUES('$sname', '$email','$organization',
	                        	   '$platform','$version','$country','$mailoption','$direction');";
        	$res = pg_query($sql_query);

//      echo "sname	   $sname<br>";       
//	echo "email	 $email<br>";	
//	echo "organization $organization<br>";
//	echo "platform     $platform<br>";	  
//	echo "version	 $version<br>"; 
//	echo "country	 $country<br>"; 
        	if ( FALSE === $res ) {
                	print "Data entry error: " . pg_last_error() . "<br />";
                	print "Contact xmipp@cnb.csic.es". "<br />";
                	print $sql_query . "<br />";
        	} else {
		print "<div align=\"center\" class=\"bigmaincell\"> Click 
		<a href=\"Xmipp-" . $version . "-" . $platform . ".tar.gz\">here</a>"; 

?>
  to dowload Xmipp and <a
  href="http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/InstallingTheSoftware#Installing_Xmipp">
  here</a> for installation instructions</div>
<?php
	        }
	}
}


if(!isset($_POST['send']) || $errmsg != '')
{
?>
<div align="center" class="errmsg"><?=$errmsg;?>
<H1>Download Xmipp</H1>
</div>
<form  method="post" name="msgform" id="msgform">
  <table width="500" border="0" align="center" cellpadding="2" cellspacing="1" class="maincell">


 <tr><td colspan="2">
<div class="note">
Please help the development of Xmipp by filling out this form correctly.
We will only use this information to justify our funding, i.e. we will not pass
it to any third parties other than our funding agencies.  
</div></td></tr>
   <tr> 
      <td width="106">Your Full Name</td>
      <td width="381"><input name="sname" type="text" class="box" id="sname" size="30" ></td>
    </tr>
    <tr> 
    <tr> 
      <td>Your Email</td>
      <td><input name="email" type="text" class="box" id="email" size="30""> <br>
    </tr>
    <tr> 
      <td>
      </td>
      <td>
      <input type="radio" name="mailoption" value="nomail"> Dont send me any mail<br>
      <input type="radio" name="mailoption" value="mail" checked> Tell me about future releases<br>
      </td>
    </tr>

      <td>Your Organization</td>
      <td><input name="organization" type="text" class="box" id="organization" size="30""></td>
    </tr>
    <tr> 
      <td>Country</td>
      <td><select name="country">

<option value="" selected="selected">Select a Country</option>

<option value="Afghanistan">Afghanistan</option>

<option value="Albania">Albania</option>

<option value="Algeria">Algeria</option>

<option value="American Samoa">American Samoa</option>

<option value="Andorra">Andorra</option>

<option value="Angola">Angola</option>

<option value="Anguilla">Anguilla</option>

<option value="Antarctica">Antarctica</option>

<option value="Antigua and Barbuda">Antigua and Barbuda</option>

<option value="Argentina">Argentina</option>

<option value="Armenia">Armenia</option>

<option value="Aruba">Aruba</option>

<option value="Australia">Australia</option>

<option value="Austria">Austria</option>

<option value="Azerbaijan">Azerbaijan</option>

<option value="Bahamas">Bahamas</option>

<option value="Bahrain">Bahrain</option>

<option value="Bangladesh">Bangladesh</option>

<option value="Barbados">Barbados</option>

<option value="Belarus">Belarus</option>

<option value="Belgium">Belgium</option>

<option value="Belize">Belize</option>

<option value="Benin">Benin</option>

<option value="Bermuda">Bermuda</option>

<option value="Bhutan">Bhutan</option>

<option value="Bolivia">Bolivia</option>

<option value="Bosnia and Herzegovina">Bosnia and Herzegovina</option>

<option value="Botswana">Botswana</option>

<option value="Bouvet Island">Bouvet Island</option>

<option value="Brazil">Brazil</option>

<option value="British Indian Ocean Territory">British Indian Ocean Territory</option>

<option value="Brunei Darussalam">Brunei Darussalam</option>

<option value="Bulgaria">Bulgaria</option>

<option value="Burkina Faso">Burkina Faso</option>

<option value="Burundi">Burundi</option>

<option value="Cambodia">Cambodia</option>

<option value="Cameroon">Cameroon</option>

<option value="Canada">Canada</option>

<option value="Cape Verde">Cape Verde</option>

<option value="Cayman Islands">Cayman Islands</option>

<option value="Central African Republic">Central African Republic</option>

<option value="Chad">Chad</option>

<option value="Chile">Chile</option>

<option value="China">China</option>

<option value="Christmas Island">Christmas Island</option>

<option value="Cocos (Keeling) Islands">Cocos (Keeling) Islands</option>

<option value="Colombia">Colombia</option>

<option value="Comoros">Comoros</option>

<option value="Congo">Congo</option>

<option value="Congo, The Democratic Republic of The">Congo, The Democratic Republic of The</option>

<option value="Cook Islands">Cook Islands</option>

<option value="Costa Rica">Costa Rica</option>

<option value="Cote D'ivoire">Cote D'ivoire</option>

<option value="Croatia">Croatia</option>

<option value="Cuba">Cuba</option>

<option value="Cyprus">Cyprus</option>

<option value="Czech Republic">Czech Republic</option>

<option value="Denmark">Denmark</option>

<option value="Djibouti">Djibouti</option>

<option value="Dominica">Dominica</option>

<option value="Dominican Republic">Dominican Republic</option>

<option value="Ecuador">Ecuador</option>

<option value="Egypt">Egypt</option>

<option value="El Salvador">El Salvador</option>

<option value="Equatorial Guinea">Equatorial Guinea</option>

<option value="Eritrea">Eritrea</option>

<option value="Estonia">Estonia</option>

<option value="Ethiopia">Ethiopia</option>

<option value="Falkland Islands (Malvinas)">Falkland Islands (Malvinas)</option>

<option value="Faroe Islands">Faroe Islands</option>

<option value="Fiji">Fiji</option>

<option value="Finland">Finland</option>

<option value="France">France</option>

<option value="French Guiana">French Guiana</option>

<option value="French Polynesia">French Polynesia</option>

<option value="French Southern Territories">French Southern Territories</option>

<option value="Gabon">Gabon</option>

<option value="Gambia">Gambia</option>

<option value="Georgia">Georgia</option>

<option value="Germany">Germany</option>

<option value="Ghana">Ghana</option>

<option value="Gibraltar">Gibraltar</option>

<option value="Greece">Greece</option>

<option value="Greenland">Greenland</option>

<option value="Grenada">Grenada</option>

<option value="Guadeloupe">Guadeloupe</option>

<option value="Guam">Guam</option>

<option value="Guatemala">Guatemala</option>

<option value="Guinea">Guinea</option>

<option value="Guinea-bissau">Guinea-bissau</option>

<option value="Guyana">Guyana</option>

<option value="Haiti">Haiti</option>

<option value="Heard Island and Mcdonald Islands">Heard Island and Mcdonald Islands</option>

<option value="Holy See (Vatican City State)">Holy See (Vatican City State)</option>

<option value="Honduras">Honduras</option>

<option value="Hong Kong">Hong Kong</option>

<option value="Hungary">Hungary</option>

<option value="Iceland">Iceland</option>

<option value="India">India</option>

<option value="Indonesia">Indonesia</option>

<option value="Iran, Islamic Republic of">Iran, Islamic Republic of</option>

<option value="Iraq">Iraq</option>

<option value="Ireland">Ireland</option>

<option value="Israel">Israel</option>

<option value="Italy">Italy</option>

<option value="Jamaica">Jamaica</option>

<option value="Japan">Japan</option>

<option value="Jordan">Jordan</option>

<option value="Kazakhstan">Kazakhstan</option>

<option value="Kenya">Kenya</option>

<option value="Kiribati">Kiribati</option>

<option value="Korea, Democratic People's Republic of">Korea, Democratic People's Republic of</option>

<option value="Korea, Republic of">Korea, Republic of</option>

<option value="Kuwait">Kuwait</option>

<option value="Kyrgyzstan">Kyrgyzstan</option>

<option value="Lao People's Democratic Republic">Lao People's Democratic Republic</option>

<option value="Latvia">Latvia</option>

<option value="Lebanon">Lebanon</option>

<option value="Lesotho">Lesotho</option>

<option value="Liberia">Liberia</option>

<option value="Libyan Arab Jamahiriya">Libyan Arab Jamahiriya</option>

<option value="Liechtenstein">Liechtenstein</option>

<option value="Lithuania">Lithuania</option>

<option value="Luxembourg">Luxembourg</option>

<option value="Macao">Macao</option>

<option value="Macedonia, The Former Yugoslav Republic of">Macedonia, The Former Yugoslav Republic of</option>

<option value="Madagascar">Madagascar</option>

<option value="Malawi">Malawi</option>

<option value="Malaysia">Malaysia</option>

<option value="Maldives">Maldives</option>

<option value="Mali">Mali</option>

<option value="Malta">Malta</option>

<option value="Marshall Islands">Marshall Islands</option>

<option value="Martinique">Martinique</option>

<option value="Mauritania">Mauritania</option>

<option value="Mauritius">Mauritius</option>

<option value="Mayotte">Mayotte</option>

<option value="Mexico">Mexico</option>

<option value="Micronesia, Federated States of">Micronesia, Federated States of</option>

<option value="Moldova, Republic of">Moldova, Republic of</option>

<option value="Monaco">Monaco</option>

<option value="Mongolia">Mongolia</option>

<option value="Montserrat">Montserrat</option>

<option value="Morocco">Morocco</option>

<option value="Mozambique">Mozambique</option>

<option value="Myanmar">Myanmar</option>

<option value="Namibia">Namibia</option>

<option value="Nauru">Nauru</option>

<option value="Nepal">Nepal</option>

<option value="Netherlands">Netherlands</option>

<option value="Netherlands Antilles">Netherlands Antilles</option>

<option value="New Caledonia">New Caledonia</option>

<option value="New Zealand">New Zealand</option>

<option value="Nicaragua">Nicaragua</option>

<option value="Niger">Niger</option>

<option value="Nigeria">Nigeria</option>

<option value="Niue">Niue</option>

<option value="Norfolk Island">Norfolk Island</option>

<option value="Northern Mariana Islands">Northern Mariana Islands</option>

<option value="Norway">Norway</option>

<option value="Oman">Oman</option>

<option value="Pakistan">Pakistan</option>

<option value="Palau">Palau</option>

<option value="Palestinian Territory, Occupied">Palestinian Territory, Occupied</option>

<option value="Panama">Panama</option>

<option value="Papua New Guinea">Papua New Guinea</option>

<option value="Paraguay">Paraguay</option>

<option value="Peru">Peru</option>

<option value="Philippines">Philippines</option>

<option value="Pitcairn">Pitcairn</option>

<option value="Poland">Poland</option>

<option value="Portugal">Portugal</option>

<option value="Puerto Rico">Puerto Rico</option>

<option value="Qatar">Qatar</option>

<option value="Reunion">Reunion</option>

<option value="Romania">Romania</option>

<option value="Russian Federation">Russian Federation</option>

<option value="Rwanda">Rwanda</option>

<option value="Saint Helena">Saint Helena</option>

<option value="Saint Kitts and Nevis">Saint Kitts and Nevis</option>

<option value="Saint Lucia">Saint Lucia</option>

<option value="Saint Pierre and Miquelon">Saint Pierre and Miquelon</option>

<option value="Saint Vincent and The Grenadines">Saint Vincent and The Grenadines</option>

<option value="Samoa">Samoa</option>

<option value="San Marino">San Marino</option>

<option value="Sao Tome and Principe">Sao Tome and Principe</option>

<option value="Saudi Arabia">Saudi Arabia</option>

<option value="Senegal">Senegal</option>

<option value="Serbia and Montenegro">Serbia and Montenegro</option>

<option value="Seychelles">Seychelles</option>

<option value="Sierra Leone">Sierra Leone</option>

<option value="Singapore">Singapore</option>

<option value="Slovakia">Slovakia</option>

<option value="Slovenia">Slovenia</option>

<option value="Solomon Islands">Solomon Islands</option>

<option value="Somalia">Somalia</option>

<option value="South Africa">South Africa</option>

<option value="South Georgia and The South Sandwich Islands">South Georgia and The South Sandwich Islands</option>

<option value="Spain">Spain</option>

<option value="Sri Lanka">Sri Lanka</option>

<option value="Sudan">Sudan</option>

<option value="Suriname">Suriname</option>

<option value="Svalbard and Jan Mayen">Svalbard and Jan Mayen</option>

<option value="Swaziland">Swaziland</option>

<option value="Sweden">Sweden</option>

<option value="Switzerland">Switzerland</option>

<option value="Syrian Arab Republic">Syrian Arab Republic</option>

<option value="Taiwan, Province of China">Taiwan, Province of China</option>

<option value="Tajikistan">Tajikistan</option>

<option value="Tanzania, United Republic of">Tanzania, United Republic of</option>

<option value="Thailand">Thailand</option>

<option value="Timor-leste">Timor-leste</option>

<option value="Togo">Togo</option>

<option value="Tokelau">Tokelau</option>

<option value="Tonga">Tonga</option>

<option value="Trinidad and Tobago">Trinidad and Tobago</option>

<option value="Tunisia">Tunisia</option>

<option value="Turkey">Turkey</option>

<option value="Turkmenistan">Turkmenistan</option>

<option value="Turks and Caicos Islands">Turks and Caicos Islands</option>

<option value="Tuvalu">Tuvalu</option>

<option value="Uganda">Uganda</option>

<option value="Ukraine">Ukraine</option>

<option value="United Arab Emirates">United Arab Emirates</option>

<option value="United Kingdom">United Kingdom</option>

<option value="United States">United States</option>

<option value="United States Minor Outlying Islands">United States Minor Outlying Islands</option>

<option value="Uruguay">Uruguay</option>

<option value="Uzbekistan">Uzbekistan</option>

<option value="Vanuatu">Vanuatu</option>

<option value="Venezuela">Venezuela</option>

<option value="Viet Nam">Viet Nam</option>

<option value="Virgin Islands, British">Virgin Islands, British</option>

<option value="Virgin Islands, U.S.">Virgin Islands, U.S.</option>

<option value="Wallis and Futuna">Wallis and Futuna</option>

<option value="Western Sahara">Western Sahara</option>

<option value="Yemen">Yemen</option>

<option value="Zambia">Zambia</option>

<option value="Zimbabwe">Zimbabwe</option>

</select></td>
    </tr>
<tr> 
  <td>Xmipp Version</td>
  <td><select name="version">
       <option value="" selected="selected">Select Xmipp Version</option>
       <option value="2.0.2">Xmipp 2.0.2</option>
       <option value="2.2">Xmipp 2.2</option>
       <option value="2.3">Xmipp 2.3</option>
       <option value="2.4">Xmipp 2.4</option>
  </select></td>
</tr>
<tr> 
  <td>Platform</td>
  <td><select name="platform">
       <option value="" selected="selected">Choose Platform</option>
       <option value="i386">linux 32 bits intel</option>
       <option value="x64">linux 64 bits intel</option>
       <option value="src">source code</option>
 </select></td>
</tr>
<tr><td colspan="2">
</td></tr>
    <tr align="center"> 
      <td colspan="2"><input name="send" type="submit" class="bluebox" id="send" value="Download" onclick="return checkForm();"></td>
    </tr>
    <tr align="center"> 
      <td colspan="2">&nbsp;</td>
    </tr>	
    <tr align="left"> 
      <td colspan="2">If by any chance the form isn't working you may contact 
        us on 
		<script language="JavaScript">
		var addr = 'xmipp';
		var host = 'cnb.csic.es';
		var email = '<a href="mailto:' + addr + '@' + host + '">' + addr + '@' + host + '</a>';
		document.write(email);
		</script></td>
    </tr>
  </table>
</form>
<?php
}

function isEmail($email)
{
	return(preg_match("/^[-_.[:alnum:]]+@((([[:alnum:]]|[[:alnum:]][[:alnum:]-]*[[:alnum:]])\.)+(ad|ae|aero|af|ag|ai|al|am|an|ao|aq|ar|arpa|as|at|au|aw|az|ba|bb|bd|be|bf|bg|bh|bi|biz|bj|bm|bn|bo|br|bs|bt|bv|bw|by|bz|ca|cc|cd|cf|cg|ch|ci|ck|cl|cm|cn|co|com|coop|cr|cs|cu|cv|cx|cy|cz|de|dj|dk|dm|do|dz|ec|edu|ee|eg|eh|er|es|et|eu|fi|fj|fk|fm|fo|fr|ga|gb|gd|ge|gf|gh|gi|gl|gm|gn|gov|gp|gq|gr|gs|gt|gu|gw|gy|hk|hm|hn|hr|ht|hu|id|ie|il|in|info|int|io|iq|ir|is|it|jm|jo|jp|ke|kg|kh|ki|km|kn|kp|kr|kw|ky|kz|la|lb|lc|li|lk|lr|ls|lt|lu|lv|ly|ma|mc|md|mg|mh|mil|mk|ml|mm|mn|mo|mp|mq|mr|ms|mt|mu|museum|mv|mw|mx|my|mz|na|name|nc|ne|net|nf|ng|ni|nl|no|np|nr|nt|nu|nz|om|org|pa|pe|pf|pg|ph|pk|pl|pm|pn|pr|pro|ps|pt|pw|py|qa|re|ro|ru|rw|sa|sb|sc|sd|se|sg|sh|si|sj|sk|sl|sm|sn|so|sr|st|su|sv|sy|sz|tc|td|tf|tg|th|tj|tk|tm|tn|to|tp|tr|tt|tv|tw|tz|ua|ug|uk|um|us|uy|uz|va|vc|ve|vg|vi|vn|vu|wf|ws|ye|yt|yu|za|zm|zw)$|(([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5])\.){3}([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5]))$/i"
			,$email));
}
?>

</body>
</html>
