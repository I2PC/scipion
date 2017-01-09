 /*****************************************************************************
 *
 * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'scipion@cnb.csic.es'
 *
 ******************************************************************************/

$(document).ready(function() {
	/*	
	* Overray the post simple method in the protocol form template. 
	* Depend a variable of the protocol form 
	*/
	$("#downloadForm").submit(function() {
		
		fullName = $("input[name=fullName]").val()
		organization = $("input[name=organization]").val()
		email = $("input[name=email]").val()
		country = $("select[name=country]").val()
		version = $("select[name=version]").val()
		platform = $("select[name=platform]").val()
		
		errors = ""
		if (fullName.length == 0)
			errors += "Please fill in the Full Name field.<br />"
		if (organization.length == 0)
			errors += "Please fill in the Organization field.<br />"
		if (email.length == 0)
			errors += "Please fill in the Email field.<br />"
		if (country.length == 0)
			errors += "Please fill in the Country field.<br />"
		if (version.length == 0)
			errors += "Please fill in the Version field.<br />"
		if (platform.length == 0)
			errors += "Please fill in the Platform field.<br />"
        if (errors.length > 0) {
        		// Show errors in the validation
        		errorPopup('Errors found', errors);
        		// Important. Stop the normal POST
        		return false
        } 
        else
        {
			var action = "/doDownload/";
			var URL = getSubDomainURL() + action
			
			var serialize_form = $(this).serialize();
			
			$.post(URL, serialize_form, function(json) {
			}, "json");
        }	
		
	});
});	