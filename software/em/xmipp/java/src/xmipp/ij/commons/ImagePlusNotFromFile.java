/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

import ij.ImagePlus;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;

/**
 *
 * @author airen
 */
public class ImagePlusNotFromFile extends ImagePlusReader{


    public ImagePlusNotFromFile(ImagePlus imp, ImageGeneric ig) {
        if(imp == null && ig == null)
            throw new IllegalArgumentException(XmippMessage.getEmptyFieldMsg("image") );
        this.imp = imp;
        this.ig = ig;
    }

  

    @Override
    public boolean getAllowsPoll() {
        return false;
    }



    @Override
    public String getName() {

        if(imp == null && ig == null)
            return null;
        String name = (ig != null)? Filename.getBaseName(ig.getFilename()): imp.getTitle();
        
        if(index != -2)
            name = String.format("%d@%s", index, name);
        return name;
    }


}
