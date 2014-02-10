/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

import ij.ImagePlus;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;

/**
 *
 * @author airen
 */
public class ImagePlusNotFromFile extends ImagePlusReader{

    ImagePlusNotFromFile(ImagePlus imp, ImageGeneric ig) {
        this.imp = imp;
        this.ig = ig;
    }

  

    @Override
    public boolean getAllowsPoll() {
        return false;
    }



    @Override
    public String getName() {
        if(imp != null)
            return imp.getTitle();
        return null;
    }
}
