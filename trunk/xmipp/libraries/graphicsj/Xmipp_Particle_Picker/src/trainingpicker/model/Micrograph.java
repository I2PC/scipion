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
		this.name = getName(file);
		this.outputfilename = name + ext;
		
	}
	
	public static String getName(String file)
	{
		String[] tokens = file.split(File.separator);
		if(tokens.length < 2)
			throw new IllegalArgumentException("Name for micrograph" +
					"is taken from parent dir, invalid path " + file);
		return  tokens[tokens.length - 2];
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
