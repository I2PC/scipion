/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

/** Implement some string utilities */
public class XmippStringUtils {

	/** Find the greatest common prefix of a group of strings */
	public static String commonPrefix(String[] strings){
		String prefix = "";
		
		int n = strings.length;
		boolean common = true;
		
		int minLenght = Integer.MAX_VALUE;

		for (int i = 0; i < n; ++i)
			minLenght = Math.min(minLenght, strings[i].length());
		
		for (int j = 0; j < minLenght && common; ++j){
			char c = strings[0].charAt(j);
			for (int i = 1; i < strings.length && common; ++i)
				if (c != strings[i].charAt(j))
					common = false;
			prefix += c;
		}

		return prefix;
	}//function commonPrefix
	
	/** Common prefix until the last "/" 
	 * This will only work for unix path separator "/" */
	public static String commonPathPrefix(String [] strings){
		String prefix = commonPrefix(strings);
		
		int pos = prefix.lastIndexOf('/');
		if (pos > 0) {
			++pos; //If found /, also include in prefix
			prefix = prefix.substring(0, pos);
		}
		else 
			prefix = "";
		
		return prefix;
	}
        
        public static String getFileExtension(String filename) {

            try {
                return filename.substring(filename.lastIndexOf("."));

            } catch (Exception e) {
                return "";
        }

}

}
