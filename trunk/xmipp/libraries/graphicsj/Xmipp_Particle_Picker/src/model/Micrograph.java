package model;

import ij.ImagePlus;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import xmipp.MDLabel;
import xmipp.MetaData;

public class Micrograph {
	
	private String filename;
	private String name;
	private List<Particle> particles;
	private ImagePlus image;
	private String ofilename;
	private static String ext = ".pos";
	
	public Micrograph(String filename, String name) {
		this.filename = filename;
		this.name = getName(filename);
		particles = new ArrayList<Particle>();
		this.ofilename = PPConfiguration.getOutputPath(name + ext);
	}
	
	public static String getName(String filename)
	{
		String[] tokens = filename.split(File.separator);
		if(tokens.length < 2)
			throw new IllegalArgumentException("Name for micrograph" +
					"is taken from parent dir, invalid path " + filename);
		return  tokens[tokens.length - 2];
	}
	
	public static String getIFilename()
	{
		return PPConfiguration.getMicrographsSelFile();
	}
	
	public String getOFilename()
	{
		return ofilename;
	}
	
	public ImagePlus getImage()
	{
		if(image == null)
			image = new ImagePlus(filename);
		return image;
	}
	
	public String getFilename() {
		return filename;
	}
	public void setFilename(String filename) {
		this.filename = filename;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	
	public List<Particle> getParticles() {
		return particles;
	}

	public void addParticle(Particle p)
	{
		particles.add(p);
	}
	

	public void removeParticle(Particle p)
	{
		particles.remove(p);
	}
	
	public String toString()
	{
		return name;
	}

	

}
