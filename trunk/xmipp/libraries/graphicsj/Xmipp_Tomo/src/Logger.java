/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
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

import java.util.Calendar;


// TODO: replace with Logger / ConsoleHandler
public class Logger {

	/** print s to stderr with a timestamp
	 * @param s
	 */
	public static void debug(String s){
		Calendar calendar = Calendar.getInstance();
		java.util.Date now = calendar.getTime();
		String output="" + now.getTime() + " > " + s;
		if(Xmipp_Tomo.TESTING == 1){
			/*if(tw != null)
				IJ.write(output);
			else */
				System.err.println(output);
		}
	}

	/** same as Debug, plus ex stack trace
	 * @param s
	 * @param ex Exception
	 */
	public static void debug(String s, Exception ex){
		debug(s + ex.toString());
		ex.printStackTrace();
	}

}
