/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.ij.commons;

import ij.ImagePlus;
import java.util.logging.Level;
import java.util.logging.Logger;
import xmipp.jni.ImageGeneric;

/**
 *
 * @author airen
 */
public abstract class ImagePlusReader {

    protected ImagePlus imp;
    protected ImageGeneric ig;
    protected boolean allowsPoll;

    
    protected long modified;
    
    protected long index = -1;
    protected boolean normalize;
    protected double normalize_min;
    protected double normalize_max;

    public abstract boolean getAllowsPoll();

    public ImagePlus getImagePlus() {
        if (imp == null) {
            imp = loadImagePlus();
        }
        return imp;
    }

    public ImagePlus loadImagePlus()
    {
		imp = null;
		try
		{
			if(ig != null)
                        {
                            
                             if(index == -1)
                                imp = XmippImageConverter.convertToImagePlus(ig);
                            else 
                             {
                                 if(ig.isStack())
                                    imp = XmippImageConverter.convertToImagePlus(ig, ImageGeneric.FIRST_IMAGE, (int)index);
                                else
                                    imp = XmippImageConverter.convertToImagePlus(ig, ImageGeneric.FIRST_IMAGE, (int)index);//read slice
                             }
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
    
    public void setNormalize(double normalize_min, double normalize_max)
    {
            this.normalize = true;
            this.normalize_min = normalize_min;
            this.normalize_max = normalize_max;

    }
    
    public boolean isVolume()
	{
		try
		{
                    
                    if (ig != null)
                        return ig.isVolume();
                    if(imp != null)
                        return imp.getStackSize() > 1;
                    loadImagePlus();
                    return isVolume();

		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e.getMessage());
		}
	}

    void setIndex(int index) {
        this.index = index;
    }

    public abstract String getName() ;

    boolean isStackOrVolume() {
        try {
            if(ig != null)
                return ig.isStackOrVolume();
            if (imp != null)
                return imp.getStackSize() > 1;
            return false;
        } catch (Exception ex) {
            Logger.getLogger(ImagePlusReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        return false;
    }
    
   

}
