package particlepicker;

import ij.ImagePlus;

import java.io.File;
import java.util.logging.Level;

public abstract class Micrograph
{

	private String file;
	private String name;
	private ImagePlus image;
	private String outputfilename;
	public static final String ext = ".pos";

	public Micrograph(String file)
	{
		this.file = file;
		if (!new File(file).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("file", file));
		this.name = getName(file, 1);
		
		this.outputfilename = name + ext;

	}

	public Micrograph(String file, String name)
	{
		this.file = file;
		if (!new File(file).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("file", file));
		this.name = name;
		this.outputfilename = name + ext;

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
		return outputfilename;
	}
	

	public ImagePlus getImagePlus()
	{
		try
		{

			if (image == null)
			{
				System.out.println("creating imageplus");
				if (file.endsWith(".tif"))
				{
					xmipp.io.readers.ImageReader reader = new xmipp.io.readers.ImageReader();
					reader.run(file);
					image = reader;

				}
				else
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
