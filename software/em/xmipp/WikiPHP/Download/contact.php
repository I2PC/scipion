<?php

$errmsg       = ''; // error message
$sname        = ''; // sender's name
$email        = ''; // sender's email address
$organization = ''; //
$platform     = ''; // binaries compiled for... platform
$architecture = ''; // architecture of the destination OS
$version      = ''; // xmipp version
$country      = ''; // country
$mailoption   = ''; // mailoption

if(isset($_POST['send'])){
  $sname         = $_POST['sname'];
  $email         = $_POST['email'];
  $organization  = $_POST['organization'];
  $platform      = $_POST['platform'];
  $version       = $_POST['version'];
  $country       = $_POST['country'];
  $mailoption    = $_POST['mailoption'];
  $architecture  = $_POST['architecture'];
  $direction     = $_SERVER['REMOTE_ADDR'];

  print "$sname<br/>";
  print "$email<br/>";
  print "$organization<br/>";
  print "$platform<br/>";
  print "$version<br/>";
  print "$country<br/>";
  print "$mailoption<br/>";
  print "$architecture<br/>";
	
  $conn = pg_connect("host=xxxx port=xxxx dbname=xxxx user=xxxxxxn password=xxxxxxx");

  if(trim($sname) == ''){
    $errmsg = 'Please enter your name';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(trim($organization) == '') {
    $errmsg = 'Please enter your organization';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(trim($email) == ''){
    $errmsg = 'Please enter your email address';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(!isEmail($email)){
    $errmsg = 'Your email address is not valid';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(trim($platform) == ''){
    $errmsg = 'Please select platform';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(trim($country) == ''){
    $errmsg = 'Please select your country';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(trim($version) == ''){
    $errmsg = 'Please select Xmipp version';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(trim($mailoption) == ''){
    $errmsg = 'Please select an email option';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if((trim($architecture) === '') && (trim($platform) !== 'src')){
    $errmsg = 'Please select the architecture option';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  } else if(((trim($architecture) !== 'i386') || (trim($version) !== '3.0')) && (trim($platform) === 'macosxbin') ){
    $errmsg = 'MAC OS X Xmipp is only available for i386 platform and v3.0';
    header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/NewDownload?errormsg=' . "$errmsg" ) ;
  }

  if($errmsg == ''){
    if(get_magic_quotes_gpc()){
      $sname        = stripslashes($sname);
      $email        = stripslashes($email);
      $organization = stripslashes($organization);
      $platform     = stripslashes($platform);
      $version      = stripslashes($version);
      $country      = stripslashes($country);
      $mailoption   = stripslashes($mailoption);
      $architecture = stripslashes($architecture);
      $direction    = stripslashes($direction);
    }	
    $sname        = utf8_encode(pg_escape_string($sname));
    $email        = pg_escape_string($email);
    $organization = utf8_encode(pg_escape_string($organization));
    $platform     = pg_escape_string($platform);
    $version      = pg_escape_string($version);
    $country      = utf8_encode(pg_escape_string($country));
    $mailoption   = pg_escape_string($mailoption);
    $architecture = pg_escape_string($architecture);
    $direction    = pg_escape_string($direction);

    // extra validations
    // the source code works for every architecture
    if($platform==='src') $architecture='all';

    //Data insertion into the xmipp_users table (xmipp database)
    $sql_query = "INSERT INTO xmipp_users (sname,email,organization,platform,version,country,mailoption,direction,os) 
    VALUES('$sname', '$email','$organization','$architecture','$version','$country','$mailoption','$direction','$platform');";
    $res = pg_query($sql_query);

    if ( FALSE === $res ){
      print "Data entry error: " . pg_last_error() . "<br />";
      print "Contact xmipp@cnb.uam.es". "<br />";
      print $sql_query . "<br />";
    } else{
      if ($platform != 'src') {
        if($platform == 'linuxbin' && $architecture == 'i386'){
          switch($version){
            case "3.0":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp3Linuxbini386Download' ) ;
              break;
            case "2.4":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp24Linuxbini386Download' ) ;
              break;
            case "2.3":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp23Linuxbini386Download' ) ;
              break;
            case "2.2":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp22Linuxbini386Download' ) ;
              break;
            case "2.0.2":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp202Linuxbini386Download' ) ;
              break;
          }
        } else if($platform == 'linuxbin' && $architecture == 'x64'){
          switch($version){
            case "3.0":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp3Linuxbinx64Download' ) ;
              break;
            case "2.4":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp24Linuxbinx64Download' ) ;
              break;
            case "2.3":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp23Linuxbinx64Download' ) ;
              break;
            case "2.2":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp22Linuxbinx64Download' ) ;
              break;
            case "2.0.2":
              header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp202Linuxbinx64Download' ) ;
              break;
          }
        } else if($platform == 'macosxbin'){
            header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp3Macosxbini386Download' ) ;
        }
        print "<div align=\"center\" class=\"bigmaincell\"> Click <a href=\"Xmipp-" . $version . "-" . $platform . "-" . $architecture . ".tar.gz\">here</a>"; 
      } else{
        switch($version){
          case "3.0":
            header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp30SrcDownload' ) ;
            break;
          case "2.4":
            header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp24SrcDownload' ) ;
            break;
          case "2.3":
            header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp23SrcDownload' ) ;
            break;
          case "2.2":
            header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp22SrcDownload' ) ;
            break;
          case "2.0.2":
            header( 'Location: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Xmipp202SrcDownload' ) ;
            break;
        }
        print "<div align=\"center\" class=\"bigmaincell\"> Click <a href=\"Xmipp-" . $version . "-" . $platform . ".tar.gz\">here</a>"; 
      }
?>
  to dowload Xmipp and <a href="http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/HowToInstall">
  here</a> for installation instructions</div>
<?php
    }
  } 
}


if(!isset($_POST['send']) || $errmsg != ''){
?>
<div align="center" class="errmsg"><?=$errmsg;?>
<!-- <H1>Download Xmipp</H1> -->
</div>
<form  action="contact.php" method="post" name="msgform" id="msgform">
  <table width="500" border="0" align="center" cellpadding="2" cellspacing="1" class="twikiTable">
    <tr class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0"> 
      <td width="306" bgcolor="#ffffff" valign="top" class="twikiTableCol0 twikiFirstCol">Your Full Name</td>
      <td width="381" bgcolor="#ffffff" valign="top" class="twikiTableCol1 twikiLastCol"><input name="sname" type="text" class="box" id="sname" size="30" ></td>
    </tr>
    <tr class="twikiTableEven twikiTableRowdataBgSorted1 twikiTableRowdataBg1">
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol0 twikiFirstCol">Your Organization</td>
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol1 twikiLastCol"><input name="organization" type="text" class="box" id="organization" size="30""></td>
    </tr>
    <tr class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0"> 
      <td bgcolor="#ffffff" valign="top" class="twikiTableCol0 twikiFirstCol">Your Email</td>
      <td bgcolor="#ffffff" valign="top" class="twikiTableCol1 twikiLastCol"><input name="email" type="text" class="box" id="email" size="30""> <br>
    </tr>
    <tr class="twikiTableEven twikiTableRowdataBgSorted1 twikiTableRowdataBg1"> 
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol0 twikiFirstCol"></td>
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol1 twikiLastCol">
      <input type="radio" name="mailoption" value="nomail"> Don't subscribe me to Xmipp users' mail list<br>
      <input type="radio" name="mailoption" value="mail" checked> Subscribe me to Xmipp users' mail list<br>
      </td>
    </tr>
    <tr class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0"> 
      <td bgcolor="#ffffff" valign="top" class="twikiTableCol0 twikiFirstCol">Country</td>
      <td bgcolor="#ffffff" valign="top" class="twikiTableCol1 twikiLastCol">
        <select name="country">
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
        </select>
      </td>
    </tr>
    <tr class="twikiTableEven twikiTableRowdataBgSorted1 twikiTableRowdataBg1"> 
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol0 twikiFirstCol">Xmipp Version</td>
      <td bgcolor="#edf4f9" valign="top" class="twikiTableCol1 twikiLastCol">
        <select name="version">
          <option value="" selected="selected">Select Xmipp Version</option>
          <option value="3.0">Xmipp 3.0</option>
          <option value="2.4">Xmipp 2.4</option>
          <option value="2.3">Xmipp 2.3</option>
          <option value="2.2">Xmipp 2.2</option>
          <option value="2.0.2">Xmipp 2.0.2</option>
        </select>
      </td>
    </tr>
  <tr class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0"> 
    <td bgcolor="#ffffff" valign="top" class="twikiTableCol0 twikiFirstCol">Platform</td>
    <td bgcolor="#ffffff" valign="top" class="twikiTableCol1 twikiLastCol">
      <select name="platform">
        <option value="" selected="selected">Choose Platform</option>
<!--        <option value="debian">Linux Ubuntu 11.10 - 12.04/Debian 6.0.4/Mint 12</option> -->
<!--        <option value="suse">Linux OpenSuSE 12.1</option> -->
<!--        <option value="fedora">Linux Fedora 16</option> -->
<!--        <option value="arch">Linux Arch</option> -->
<!--        <option value="gentoo">Linux Gentoo</option> -->
<!--        <option value="win">Windows</option> -->
        <option value="src">Generic (Source Code)</option>
        <option value="linuxbin">Linux binaries</option>
        <option value="macosxbin">MAC OS X binaries</option>
      </select>
    </td>
  </tr>
  <tr class="twikiTableEven twikiTableRowdataBgSorted1 twikiTableRowdataBg1">
    <td bgcolor="#edf4f9" valign="top" class="twikiTableCol0 twikiFirstCol">Architecture</td>
    <td bgcolor="#edf4f9" valign="top" class="twikiTableCol1 twikiLastCol">
      <select name="architecture">
        <option value="" selected="selected">Choose Architecture</option>
        <option value="i386">32 bits</option>
        <option value="x64">64 bits</option>
      </select>
    </td>
  </tr> 
  <tr align="center" class="twikiTableOdd twikiTableRowdataBgSorted0 twikiTableRowdataBg0">
    <td colspan="2"><input name="send" type="submit" class="twikiSubmit" id="send" value="Download" onclick="return checkForm();"></td>
  </tr>
  </table>
</form>
If by any chance the form isn't working you may contact us on <a href="mailto:xmipp@cnb.csic.es"> xmipp@cnb.csic.es </a>
<?php
}

function isEmail($email){
  return(preg_match("/^[-_.[:alnum:]]+@((([[:alnum:]]|[[:alnum:]][[:alnum:]-]*[[:alnum:]])\.)+(ad|ae|aero|af|ag|ai|al|am|an|ao|aq|ar|arpa|as|at|au|aw|az|ba|bb|bd|be|bf|bg|bh|bi|biz|bj|bm|bn|bo|br|bs|bt|bv|bw|by|bz|ca|cc|cd|cf|cg|ch|ci|ck|cl|cm|cn|co|com|coop|cr|cs|cu|cv|cx|cy|cz|de|dj|dk|dm|do|dz|ec|edu|ee|eg|eh|er|es|et|eu|fi|fj|fk|fm|fo|fr|ga|gb|gd|ge|gf|gh|gi|gl|gm|gn|gov|gp|gq|gr|gs|gt|gu|gw|gy|hk|hm|hn|hr|ht|hu|id|ie|il|in|info|int|io|iq|ir|is|it|jm|jo|jp|ke|kg|kh|ki|km|kn|kp|kr|kw|ky|kz|la|lb|lc|li|lk|lr|ls|lt|lu|lv|ly|ma|mc|md|mg|mh|mil|mk|ml|mm|mn|mo|mp|mq|mr|ms|mt|mu|museum|mv|mw|mx|my|mz|na|name|nc|ne|net|nf|ng|ni|nl|no|np|nr|nt|nu|nz|om|org|pa|pe|pf|pg|ph|pk|pl|pm|pn|pr|pro|ps|pt|pw|py|qa|re|ro|ru|rw|sa|sb|sc|sd|se|sg|sh|si|sj|sk|sl|sm|sn|so|sr|st|su|sv|sy|sz|tc|td|tf|tg|th|tj|tk|tm|tn|to|tp|tr|tt|tv|tw|tz|ua|ug|uk|um|us|uy|uz|va|vc|ve|vg|vi|vn|vu|wf|ws|ye|yt|yu|za|zm|zw)$|(([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5])\.){3}([0-9][0-9]?|[0-1][0-9][0-9]|[2][0-4][0-9]|[2][5][0-5]))$/i",$email));
}
?>

