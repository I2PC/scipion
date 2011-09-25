package trainingpicker.model;

import ij.ImagePlus;

import java.io.File;
import java.util.List;

public abstract class Micrograph {
	
	private String file;
	private String name;
	private ImagePlus image;
	private String outputfilename;
	private static String ext = ".pos";
	
	public Micrograph(String file) {
		this.file = file;
		if(!new File(file).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("file", file));
		this.name = getName(file, 2);
		this.outputfilename = name + ext;
		
	}
	
	public Micrograph(String file, String name) {
		this.file = file;
		if(!new File(file).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("file", file));
		this.name = name;
		this.outputfilename = name + ext;
		
	}
	
	public static String getName(String file, int level)
	{
		String[] tokens = file.split(File.separator);
		if(tokens.length < level )
			throw new IllegalArgumentException(String.format("Name for micrograph is taken from level %s, invalid path ", level, file));
		String name = tokens[tokens.length - level];
		int pos = name.lastIndexOf('.');
		name = name.substring(0, pos);
		return  name;
	}

	public String getOFilename()
	{
		return outputfilename;
	}
	
	
	
	public ImagePlus getImage()
	{
		if(image == null)
			image = new ImagePlus(file);
		return image;
	}
	
	public void releaseImage()
	{
		image = null;
	}
	
	
	public String getFile() {
		return file;
	}
	
	public String getName() {
		return name;
	}
	
		
	
	public String toString()
	{
		return name;
	}

	
	public abstract boolean hasData();
	
	public abstract void reset();


}
