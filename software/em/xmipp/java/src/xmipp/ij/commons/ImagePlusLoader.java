/**
 * *************************************************************************
 * Authors: J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * All comments concerning this program package may be sent to the e-mail
 * address 'xmipp@cnb.csic.es'
 * *************************************************************************
 */
package xmipp.ij.commons;

import ij.ImagePlus;
import java.io.File;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;

public class ImagePlusLoader {

    protected ImagePlusReader impreader;
    private boolean existsfile;
    
    public ImagePlusLoader(ImagePlus imp) {
        this(getFile(imp), imp, null, false, false, false, -2);
    }

    public ImagePlusLoader(ImagePlus imp, boolean inverty) {
        this(getFile(imp), imp, null, false, false, inverty, -2);
    }
    
    public ImagePlusLoader(ImageGeneric ig, boolean inverty) {
        this(getFile(ig), null, ig, false, false, inverty, -2);
    }

    public ImagePlusLoader(int index, ImageGeneric ig) {
        this(getFile(ig), null, ig, false, false, false, index);
    }
    
    public ImagePlusLoader(int index, ImageGeneric ig, boolean inverty) {
        this(getFile(ig), null, ig, false, false, inverty, index);
    }

    public ImagePlusLoader(String fileName, boolean inverty) {
        this(fileName, null, null, false, false, inverty, -2);
    }
    
    public ImagePlusLoader(String fileName) {
        this(fileName, null, null, false, false, false, -2);
    }

    public ImagePlusLoader(String fileName, boolean useGeometry, boolean wrap, boolean inverty) {
        this(fileName, null, null, useGeometry, wrap, inverty, -2);
    }

    public ImagePlusLoader(String fileName, ImagePlus imp, ImageGeneric ig, boolean useGeo, boolean wrap, boolean inverty, int index) {
        
        if (fileName != null) {
            String path = Filename.findImagePath(fileName, null, true);//check if file exists dismissing preffix and suffix
            existsfile = path != null;
        }
        if (existsfile) 
            impreader = new ImagePlusFromFile(fileName, imp, ig);
        else 
        {

            if(fileName != null)
                throw new IllegalArgumentException(XmippMessage.getPathNotExistsMsg(fileName));
            impreader = new ImagePlusNotFromFile(imp, ig);
        }
//        
        impreader.setIndex(index);
        impreader.setUseGeometry(useGeo);
        setWrap(wrap);
        setInvertY(inverty);
    }

    public ImagePlusLoader(ImageGeneric ig) {
        this(-1, ig);
    }

    public void setNormalize(double normalize_min, double normalize_max) {
        impreader.setNormalize(normalize_min, normalize_max);

    }
    
    public void setInvertY(boolean invert) {
        impreader.setInvertY(invert);

    }

    public ImagePlus getImagePlus() {
        return impreader.getImagePlus();
    }

    public void setWrap(boolean value) {
        impreader.setWrap(value);

    }

    public boolean isWrap() {
        return impreader.isWrap();
    }

    public static String getFile(ImagePlus imp) {
        String file = null;

        if (imp != null && imp.getOriginalFileInfo() != null) {
            file = imp.getOriginalFileInfo().directory + File.separator + imp.getOriginalFileInfo().fileName;
        }

        return file;
    }

    public static String getFile(ImageGeneric ig) {
        String file = null;
        if (ig != null && ig.getFilename() != null) {
            file = ig.getFilename();
        }
        return file;
    }

    public boolean allowsPoll() {
        return impreader.getAllowsPoll();
    }

    public String getName() {
        return impreader.getName();
    }

    public ImagePlus loadImagePlus() {
        return impreader.loadImagePlus();
    }

    public boolean existsFile() {
        return existsfile;
    }

    public boolean isVolume() {
        return impreader.isVolume();
    }


    public boolean isStackOrVolume() {
        return impreader.isStackOrVolume();
    }
    
     public boolean getUseGeometry() {
        return impreader.getUseGeometry();
    }

     public void setUseGeometry(boolean value) {
        impreader.setUseGeometry(value);

    }
     
     public void setGeometry(Geometry geometry)
    {
        impreader.setGeometry(geometry);
    }
     
      public boolean hasGeometry()
    {
        return impreader.geometry != null;
    }

    
    public void setDimension(int width, int height)
    {
        impreader.setDimension(width, height);
    }

    public boolean isInvertY() {
        return impreader.isInvertY();
    }
    
    

}
