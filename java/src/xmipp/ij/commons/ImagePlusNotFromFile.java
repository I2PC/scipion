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
   public ImagePlus loadImagePlus()
    {
            imp = null;
            try
            {
                    if (ig != null)
                    {
                            if (index != -1)
                                    imp = XmippImageConverter.convertToImagePlus(ig, ImageGeneric.FIRST_IMAGE, index);
                            else
                                    imp = XmippImageConverter.readToImagePlus(ig);
                    }
                    if(normalize)
                    {
                            imp.getProcessor().setMinAndMax(normalize_min, normalize_max);
                            imp.updateImage();
                    }
                    return imp;
            }
            catch (Exception e)
            {
                    e.printStackTrace();
            }
            return imp;
    }

    @Override
    public String getName() {
        if(imp != null)
            return imp.getTitle();
        return null;
    }
}
