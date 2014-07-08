/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.ij.commons;

import ij.ImagePlus;
import ij.process.ImageProcessor;
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
    protected boolean useGeometry;
    
    protected long modified;
    
    protected long index = -1;
    protected boolean normalize;
    protected double normalize_min;
    protected double normalize_max;
    protected Geometry geometry;
//    protected int width, height;
    

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
                            
                        checkResizeAndGeo();
			
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
    
    public void checkResizeAndGeo()
    {
        if(useGeometry && geometry != null)
        {
            try {
                
                ImageGeneric tempig = XmippImageConverter.convertToImageGeneric(imp);
                tempig.applyGeo(geometry.shiftx, geometry.shifty, geometry.psiangle);
                imp = XmippImageConverter.convertToImagePlus(tempig);
            } catch (Exception ex) {
                Logger.getLogger(ImagePlusReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
//        if(width != 0  && height != 0)
//        {
//            ImageProcessor processor = imp.getProcessor();
//            processor.setInterpolate(true);
//            processor = processor.resize(width, height);
//            imp = new ImagePlus("", processor);
//        }
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
            
            if (imp != null)
                return imp.getStackSize() > 1;
            return false;
        } catch (Exception ex) {
            Logger.getLogger(ImagePlusReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        return false;
    }
    
   
    public void setGeometry(Geometry geometry)
    {
        this.geometry = geometry;
    }
    
    public boolean getUseGeometry() {
        return useGeometry;
    }

    public void setUseGeometry(boolean value) {
        useGeometry = value;

    }

//    void setDimension(int width, int height) {
//        this.width = width;
//        this.height = height;
//    }

}
