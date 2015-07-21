/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.ij.commons;

import ij.ImagePlus;
import java.io.File;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;

/**
 *
 * @author airen
 */
public class ImagePlusFromFile extends ImagePlusReader{
    
        
    private String fileName;
    private long modified;
    
    public ImagePlusFromFile(String fileName)
    {
        if (fileName == null || fileName.equals(""))
            throw new IllegalArgumentException("empty file");
        this.fileName = fileName;
        this.modified = new File(fileName).lastModified();
        
    }
    
    public ImagePlusFromFile(String fileName, ImagePlus imp, ImageGeneric ig)
    {
        this(fileName);
        this.imp = imp;
        this.ig = ig;
    }
        
        
    @Override
	public ImagePlus loadImagePlus()
	{
                
        imp = null;
		try
		{
			if (ig == null || hasChanged())
            {
                if(Filename.isXmippSupported(fileName))
                    try
                    {
                        ig = new ImageGeneric(fileName);//to read again file
                    }
                    catch(Exception e)
                    {
                        imp = new ImagePlus(fileName);
                    }
                else
                    imp = new ImagePlus(fileName);
                if(ig != null && !hasIndex())
                {
                    imp = XmippImageConverter.readToImagePlus(ig);
                }
                else if(ig != null)
                {
                    if(ig.isStack())
                        imp = XmippImageConverter.readToImagePlus(ig, index, ImageGeneric.ALL_SLICES);//read image or volume on indexn
                    else
                        imp = XmippImageConverter.readToImagePlus(ig, (int)index);//read image slice on volume
                }
            }
            else if(ig != null)
            {
                 if(index == -2 || !ig.isStackOrVolume())//if index is not -2 ig read it already
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
		}
		return imp;
	}
        
        
    
    public boolean hasChanged()
	{
		return new File(fileName).lastModified() > modified;
	}
        

	public String getFileName()
	{
		return fileName;
	}

    @Override
    public boolean getAllowsPoll() {
        return true;
    }

    @Override
    public String getName() {

        String name = fileName;
        if(!name.contains("@"))
        	name = Filename.getBaseName(name);
        else
        	name = Filename.getPrefix(name) + Filename.SEPARATOR + Filename.getBaseName(Filename.getFilename(name));

        if(index != -2)
            name = String.format("%d@%s", index, name);
        return name;

    }

    

        
        
        
}
