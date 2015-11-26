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
    protected boolean useGeometry;
    protected boolean wrap;
    protected long modified;
    
    protected long index = -2;
    protected boolean normalize;
    protected double normalize_min;
    protected double normalize_max;
    protected Geometry geometry;
    protected int width, height;
    protected boolean inverty;
    

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
                
                 if(!hasIndex() )
                    imp = XmippImageConverter.convertToImagePlus(ig);
                 else 
                 {
                    if(ig.isStack())
                        imp = XmippImageConverter.convertToImagePlus(ig, index);//read image or volume on index
                    else
                        imp = XmippImageConverter.convertToImagePlus(ig, ImageGeneric.FIRST_IMAGE, (int)index);//read image slice on volume
                 }
            }
            checkResizeAndGeo();
            checkInvertY();
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
                        throw new IllegalArgumentException(e);
		}
	
    }
    
    public void checkResizeAndGeo()
    {
        if(useGeometry && geometry != null)
        {
            try {
                
                ImageGeneric tempig = XmippImageConverter.convertToImageGeneric(imp);
                if (geometry.hasMatrix())
                	tempig.applyGeoMatrix(geometry.getMatrix(), wrap);
                else
                	tempig.applyGeo(geometry.shiftx, geometry.shifty, geometry.psiangle, geometry.flip, wrap, geometry.scaleFactor);
                imp = XmippImageConverter.convertToImagePlus(tempig);
            } catch (Exception ex) {
                Logger.getLogger(ImagePlusReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        ImageProcessor processor = imp.getProcessor();
        if(width != 0  && height != 0 && imp != null && processor != null)
        {
            processor.setInterpolate(true);
            processor = processor.resize(width, height);
            imp = new ImagePlus("", processor);
        }
    }
    
    
    protected void checkInvertY() {
        if(inverty)
        {
            imp.getProcessor().flipVertical();
            imp.updateImage();
        }
    }
    
    
    public void setNormalize(double normalize_min, double normalize_max)
    {
            this.normalize = true;
            this.normalize_min = normalize_min;
            this.normalize_max = normalize_max;

    }
    
    public void setInvertY(boolean inverty)
    {
            this.inverty = inverty;

    }
    
    public boolean isVolume()
    {
            try
            {
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
    
    public boolean hasIndex()
    {
        return index != -2;
    }

    public abstract String getName() ;


    boolean isStackOrVolume() {
        if(imp == null)
            loadImagePlus();
        return imp.getStackSize() > 1;
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

    public void setDimension(int width, int height) {
        this.width = width;
        this.height = height;
    }
    
    public void setWrap(boolean value) {
        wrap = value;

    }

    public boolean isWrap() {
        return wrap;
    }

    public boolean isInvertY() {
        return inverty;
    }
    

}
