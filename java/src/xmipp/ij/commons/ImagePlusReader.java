/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.ij.commons;

import ij.ImagePlus;
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
    
    protected int index = -1;
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

    public abstract ImagePlus loadImagePlus();
    
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

}
