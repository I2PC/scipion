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
 **************************************************************************
 */
package xmipp.ij.commons;

import java.io.File;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import ij.ImagePlus;

public class ImagePlusLoader {

    protected ImagePlusReader impreader;
    protected boolean wrap;
    protected boolean allowsGeometry = false;
    protected boolean useGeometry;
    private boolean existsfile;

    public ImagePlusLoader() {

    }

    public ImagePlusLoader(ImagePlus imp) {
        this(getFile(imp), imp, null);

    }

    public ImagePlusLoader(int index, ImageGeneric ig) {
        this(getFile(ig), null, ig);
        impreader.setIndex(index);

    }

    public ImagePlusLoader(String fileName)
    {
        this(fileName, null, null);
    }
    
    public ImagePlusLoader(String fileName, boolean allowsGeometry, boolean useGeometry, boolean wrap)
    {
        this(fileName, null, null);
        this.allowsGeometry = allowsGeometry;
        this.useGeometry = useGeometry;
        this.wrap = wrap;
    }
            
            
    public ImagePlusLoader(String fileName, ImagePlus imp, ImageGeneric ig) {
        if (fileName == null || fileName.equals("")) {
            throw new IllegalArgumentException("File not found");
        }
        int index = -1;
        if (fileName.contains("@")) {
            int arrobaindex = fileName.lastIndexOf("@");
            String header = fileName.substring(0, arrobaindex);
            
            int sepindex = header.lastIndexOf(File.separator);//-1 if separator does not exists
            String textindex = fileName.substring(sepindex + 1, arrobaindex);
            index = Integer.parseInt(textindex);
            fileName = fileName.substring(arrobaindex + 1);
            if(sepindex != -1)
                fileName = Filename.join(header.replace(textindex + "@", ""), fileName);
        }

        if (!new File(fileName).exists()) {
            throw new IllegalArgumentException("File not found " + fileName);
        }
        existsfile = true;
        impreader = new ImagePlusFromFile(fileName, imp, ig);
        impreader.setIndex(index);

    }

    public ImagePlusLoader(ImageGeneric ig) {
        this(-1, ig);
    }

    public void setNormalize(double normalize_min, double normalize_max) {
        impreader.setNormalize(normalize_min, normalize_max);

    }

    public ImagePlus getImagePlus() {
        return impreader.getImagePlus();
    }

    public boolean allowsGeometry() {
        return allowsGeometry;
    }

    public boolean getUseGeometry() {
        return useGeometry;
    }

    public void setUseGeometry(boolean value) {
        useGeometry = value;

    }

    public void setWrap(boolean value) {
        wrap = value;

    }

    public boolean isWrap() {
        return wrap;
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
        return impreader.allowsPoll;
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

    public void setAllowsGeometry(boolean allowsGeometry) {
        this.allowsGeometry = allowsGeometry;
    }

    public boolean isStackOrVolume() {
        return impreader.isStackOrVolume();
    }

}
