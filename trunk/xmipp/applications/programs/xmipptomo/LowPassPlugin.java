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
package xmipptomo;

import ij.gui.GenericDialog;


public class LowPassPlugin extends Plugin {
	// default value - @see xmipptomo.Plugin.collectParameters, ij.plugin.FFTFilter
	private double filterLargeDia = 40.0;
	private  double  filterSmallDia = 3.0;
	// private  int choiceIndex = 0;
	private  String choiceDia = "None";
	private  double toleranceDia = 5.0;
	private  boolean doScalingDia = true;
	private  boolean saturateDia = true;
	private  boolean displayFilter;
	
	public static String COMMAND="Bandpass Filter...";
	
	public String getCommand(){
		return COMMAND;
	}
	

	@Override
	public String getOptions(){
		return "Filter_Large Structures Down to=" + filterLargeDia +
		" Filter_Small Structures Up to=" + filterSmallDia + " Suppress Stripes:=" + choiceDia +
		" Tolerance of Direction:=" + toleranceDia + " Autoscale After Filtering=" + doScalingDia + 
		" Saturate Image when Autoscaling=" + saturateDia + " Display Filter=" + displayFilter;
	}
	
	@Override
	public void collectParameters(GenericDialog gd) {
		if(gd != null){
			filterLargeDia=gd.getNextNumber();
			filterSmallDia=gd.getNextNumber();
			// choiceIndex=(int)gd.getNextNumber();
			choiceDia=gd.getNextString();
			toleranceDia=gd.getNextNumber();
			doScalingDia=gd.getNextBoolean();
			saturateDia=gd.getNextBoolean();
			displayFilter=gd.getNextBoolean();
		}
	}

	
	public String toString(){
		return getOptions();
	}

}
