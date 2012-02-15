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

package xmipp.viewer;

import java.awt.Dimension;

import ij.ImagePlus;

/** This class will contains basic info to be used for image rendering.
 * It will contains an ImagePlus, label and some other useful info.
 */
public class ImageItem {
	
	protected ImagePlus image;
	protected String key; //This will be used as id for store in cache
	protected String label;
	//Flag to check selection
	public boolean isSelected = false;
	//Flag to mark when enabled
	public boolean isEnabled = true;
	//Flag to mark if display label
	public boolean showLabel = false;
	// Cell dimension on table
	public Dimension cellDim;
	
	//Constructor
	public ImageItem(String k, String l, ImagePlus imp){
		key = k;
		image = imp;
		label = l;
	}
	
	public ImagePlus getImage() {
		return image;
	}
	
	public String getLabel() {
		return label;
	}
	
	public String getKey() {
		return key;
	}
}
