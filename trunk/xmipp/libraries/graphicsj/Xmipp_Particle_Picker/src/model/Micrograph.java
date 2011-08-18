package model;

import ij.ImagePlus;

import java.awt.Color;
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
	private String xmd;
	
	public Micrograph(String filename, String name) {
		this.filename = filename;
		this.name = name;
		image  = new ImagePlus(filename);
		particles = new ArrayList<Particle>();
		this.xmd = name + ".xmd";
	}
	
	public String getXMD()
	{
		return xmd;
	}
	
	public ImagePlus getImage()
	{
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
