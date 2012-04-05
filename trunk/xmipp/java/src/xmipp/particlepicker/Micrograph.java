package xmipp.particlepicker;

import ij.IJ;
import ij.ImagePlus;
import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.ImageGeneric;
import xmipp.utils.XmippMessage;

import java.io.File;
import java.util.List;
import java.util.logging.Level;

public abstract class Micrograph
{

	private String file;
	private String name;
	private ImagePlus imp;	
	private String posfile;
	private String pos24file;
	public static final String ext = ".pos";
	

	
	public void setPosFileFromXmipp24(String posfile)
	{
		this.pos24file = posfile;
	}
	
	public String getPosFileFromXmipp24()
	{
		return pos24file;
	}

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
		this.posfile = name + ext;

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

	public String getPosFile()
	{
		return posfile;
	}
	

	public ImagePlus getImagePlus()
	{
		try
		{

			if (imp == null)
			{
				System.out.println("Empty Filters");
				imp = XmippImageConverter.loadImage(file);
				if(imp == null)
					imp = new ImagePlus(file);
			}
			return imp;
		}
		catch (Exception e)
		{
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	public ImagePlus getImagePlus(List<Filter> filters){
		try
		{
			if (filters.isEmpty())
				return getImagePlus();
			
			if (imp == null)
			{
				ImageGeneric ig = new ImageGeneric(file);
				ig.read(ImageGeneric.FIRST_IMAGE);
				for (Filter f: filters)
					if (f.getCommand().equals("Smooth Filter")){ //this filter should be the first applied
						ig.convert2Datatype(ImageGeneric.UChar);
						ImageGeneric igsmooth = new ImageGeneric(ImageGeneric.UChar);
						System.out.printf("width: %s height: %s\n", ig.getXDim(), ig.getYDim());
						igsmooth.resize(ig.getXDim(), ig.getYDim());
						ig.smooth(igsmooth);
						imp = XmippImageConverter.convertToImagePlus(igsmooth);
						ig.destroy();
					} 
					else 
					{
						if (imp == null) {
							imp = XmippImageConverter.convertToImagePlus(ig);
							ig.destroy();
						}
						//String macro = String.format("run(\"%s\", \"%s\");", f.getCommand(), f.getOptions());
						IJ.run(imp, f.getCommand(), f.getOptions());
					}
			}
			return imp;
		}
		catch (Exception e)
		{
			ParticlePicker.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public void releaseImage()
	{
		imp = null;
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
