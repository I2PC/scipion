/***************************************************************************
 * Authors:     Airen Zaldivar (airen@cnb.csic.es)
 *
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

package xmipp.utils;


/** This class will contains methods to return formatted text messages */
public class XmippMessage {
	

	/** Particle Picker messages */
	public static String getIllegalDeleteMsg(String item)
	{
		return String.format("There it must be at least one %s defined", item);
	}
	
	public static String getAlreadyExistsGroupNameMsg(String name)
	{
		return String.format("Group %s already exists", name);
	}

	public static String getEmptyFieldMsg(String field) {
		return String.format("Must specify %s", field);
	}

	public static String getAssociatedDataMsg(String field) {
		return String.format("%s has associated data. Can not be removed");
	}

	public static String getNoSuchFieldValueMsg(String field, Object value) {
		return String.format("No such %s %s exists", field, value);
	}

	public static String getOutOfBoundsMsg(Object o)
	{
		return String.format("%s out of bounds", o);
	}
	
	public static String getNotImplementedYetMsg()
	{
		return "Not implemented yet";
	}
	
	public static String getIllegalValueMsg(String field, Object value)
	{
		return String.format("Illegal value for %s: %s", field, value);
	}
	
	public static String getPathNotExistsMsg(String path){
		return String.format("Path '%s' doesn't exists.", path);
	}
	
	public static String getIllegalStateForOperationMsg(String field, String value){
		return String.format("Invalid operation for %s on state %s", field, value);
	}

	public static String getIllegalValueMsgWithInfo(String field, Object value, String message)
	{
		return String.format("%s=%s\n%s", field, value, message);
	}

	public static String getUnexpectedErrorMsg()
	{
		return "Unexpected error";
	}
}
