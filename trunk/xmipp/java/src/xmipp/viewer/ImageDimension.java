/***************************************************************************
 * Authors: Juanjo Vega     
 * 			J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import ij.IJ;
import xmipp.jni.ImageGeneric;

/**
 * Simple class to store image dimensions
 */
public class ImageDimension {

    private int xdim, ydim, zdim;
    private long ndim;

    public ImageDimension() {
    }

    public ImageDimension(int x){
    	this(x, x);
    }
    
    public ImageDimension(int x, int y){
    	this(x, y, 1, 1);
    }
    
    public ImageDimension(int x, int y, int z){
    	this(x, y, z, 1);
    }
    
    public ImageDimension(int x, int y, int z, long n) {
        this.xdim = x;
        this.ydim = y;
        this.zdim = z;
        this.ndim = n;
    }

    public ImageDimension(ImageGeneric image) {
        try {
            xdim = image.getXDim();
            ydim = image.getYDim();
            zdim = image.getZDim();
            ndim = image.getNDim();
        } catch (Exception ex) {
            IJ.error("Retrieving image dimensions: "+ image);
        }
    }

    public int getXDim() {
        return xdim;
    }

    public int getYDim() {
        return ydim;
    }

    public int getZDim() {
        return zdim;
    }

    public long getNDim() {
        return ndim;
    }

    public void setXDim(int x) {
        this.xdim = x;
    }

    public void setYDim(int y) {
        this.ydim = y;
    }

    public void setZDim(int z) {
        this.zdim = z;
    }

    public void setNDim(long n) {
        this.ndim = n;
    }
}
