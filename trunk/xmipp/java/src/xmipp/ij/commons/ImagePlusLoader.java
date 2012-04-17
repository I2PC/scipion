package xmipp.ij.commons;

import java.io.File;

import xmipp.jni.ImageGeneric;
import ij.ImagePlus;

public class ImagePlusLoader
{

	private String fileName = null;
	private boolean allowsPoll;
	private boolean allowsGeometry = false;
	private boolean useGeometry;
	private ImagePlus imp;
	private ImageGeneric ig;
	private long modified;

	public ImagePlusLoader(ImagePlus imp)
	{
		this.imp = imp;
	}

	public ImagePlusLoader(String fileName)
	{
		this.fileName = fileName;
		this.modified = new File(fileName).lastModified();
		allowsPoll = true;
	}

	public ImagePlusLoader(ImageGeneric ig)
	{
		this.ig = ig;
		allowsPoll = false;
	}

	public ImagePlus getImagePlus()
	{
		if (imp == null)
			imp = loadImagePlus();
		return imp;
	}

	public ImagePlus loadImagePlus()
	{
		//ImagePlus imp = null;
		try
		{
			if (fileName != null && ( hasChanged() || imp == null))
				imp = XmippImageConverter.loadImage(fileName);
			else if (ig != null)
				imp = XmippImageConverter.readToImagePlus(ig);
			return imp;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return imp;
	}
	
	public boolean hasChanged(){
		return new File(fileName).lastModified() > modified;
	}

	public String getFileName()
	{
		return fileName;
	}

	public boolean allowsPoll()
	{
		return allowsPoll;
	}

	public boolean allowsGeometry()
	{
		return allowsGeometry;
	}

	public boolean getUseGeometry()
	{
		return useGeometry;
	}

	public void setUseGeometry(boolean useGeometry)
	{
		this.useGeometry = useGeometry;
	}

}
