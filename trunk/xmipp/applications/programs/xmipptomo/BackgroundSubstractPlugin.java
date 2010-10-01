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

import ij.gui.GenericDialog;


public class BackgroundSubstractPlugin extends Plugin {
	// default values - @see xmipptomo.Plugin.collectParameters, ij.plugin.filter.BackgroundSubstracter
	private double radius = 50.0;
    private static boolean lightBackground = true;
    private static boolean createBackground;   // don't subtract background (e.g., for processing the background before subtracting)
    private static boolean useParaboloid; // use "Sliding Paraboloid" instead of rolling ball algorithm
    private static boolean doPresmooth = true; // smoothen the image before creating the background
	
	public static String COMMAND="Subtract Background...";
	
	@Override
	public String getCommand(){
		return COMMAND;
	}
	
	// adapt radius if the original image was resized 
	@Override
	public String getOptions(){
		return "Rolling Ball Radius:=" + radius + "Light Background=" + lightBackground + "Create Background (Don't Subtract)=" +
		createBackground + "Sliding Paraboloid=" + useParaboloid + "Disable Smoothing=" + doPresmooth;
	}
	
	@Override
	public void collectParameters(GenericDialog gd) {
		if(gd != null){
			radius=gd.getNextNumber();
			lightBackground=gd.getNextBoolean();
			createBackground=gd.getNextBoolean();
			useParaboloid=gd.getNextBoolean();
			doPresmooth=!gd.getNextBoolean();
		}
	}

	
	public String toString(){
		return getOptions();
	}

}
