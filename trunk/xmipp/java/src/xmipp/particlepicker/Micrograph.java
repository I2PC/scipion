package xmipp.particlepicker;

import ij.ImagePlus;
import xmipp.ij.XmippImageConverter;
import xmipp.utils.XmippMessage;

import java.io.File;
import java.util.logging.Level;

public abstract class Micrograph
{

	private String file;
	private String name;
	private ImagePlus image;
	private String ofilename;
	public static final String ext = ".pos";

	public Micrograph(String file)
	{
		this(file, getName(file, 1));
	}

	public Micrograph(String file, String name)
	{
		this.file = file;
		if (!new File(file).exists())
			throw new IllegalArgumentException(XmippMessage.getNoSuchFieldValueMsg("file", file));
		this.name = name;
		this.ofilename = name + ext;

	}

	public static String getName(String file, int level)
	{
		String[] tokens = file.split(File.separator);
		if (tokens.length < level)
			throw new IllegalArgumentException(String.format("Name for micrograph is taken from level %s, invalid path ", level, file));
		String name = tokens[tokens.length - level];
		if (level == 1)
		{
			int pos = name.lastIndexOf('.');
			if (pos != -1)
				name = name.substring(0, pos);
		}
		return name;
	}

	public String getOFilename()
	{
		return ofilename;
	}
	

	public ImagePlus getImagePlus()
	{
		try
		{

			if (image == null)
			{

				image = XmippImageConverter.loadImage(file);
				if(image == null)
					image = new ImagePlus(file);
			}
			return image;
		}
		catch (Exception e)
		{
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void releaseImage()
	{
		image = null;
	}

	public String getFile()
	{
		return file;
	}

	public String getName()
	{
		return name;
	}

	public String toString()
	{
		return name;
	}

	public abstract boolean hasData();

	public abstract void reset();

}
