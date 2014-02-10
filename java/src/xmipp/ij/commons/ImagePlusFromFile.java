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
            if (Filename.isVolume(fileName))
                try
                {
                        ig = new ImageGeneric(fileName);
                }
                catch (Exception e)
                {
                        throw new IllegalArgumentException(e.getMessage());
                }
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
                            ig = new ImageGeneric(fileName);//to read again file
                            
                        }
                        if(index == -1)
                            imp = XmippImageConverter.readToImagePlus(ig);
                        else
                        {
                             
                             imp = XmippImageConverter.readToImagePlus(ig, index);//reading 
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
        return fileName;
    }
        
        
}
