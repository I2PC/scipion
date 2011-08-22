package model;

import ij.ImagePlus;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;

public class Micrograph {
	
	private String filename;
	private String name;
	private List<Particle> particles;
	private ImagePlus image;
	private String outputfname;
	private static String ext = ".pos";
	private String ctf;
	private ImageIcon ctficon;
	private String aoutputfname;
	private boolean autopicking = false;
	
	
	public Micrograph(String filename, String ctf) {
		this.filename = filename;
		this.name = getName(filename);
		particles = new ArrayList<Particle>();
		this.outputfname = PPConfiguration.getOutputPath(name + ext);
		this.aoutputfname = PPConfiguration.getOutputPath(name + "_auto" + ext);
		this.ctf = ctf;
	}
	
	public boolean isAutopicking()
	{
		return autopicking;
	}
	
	public void setAutopicking(boolean autopicking)
	{
		this.autopicking = autopicking;
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
	
	public String getOutputFName()
	{
		return outputfname;
	}
	public String getAutoOutputFName()
	{
		return aoutputfname;
	}
	
	public ImagePlus getImage()
	{
		if(image == null)
			image = new ImagePlus(filename);
		return image;
	}
	
	public Icon getCTFIcon()
	{
		if(ctficon == null)
		{
			ImagePlus ip = new ImagePlus(ctf);
			//ip.getProcessor().scale(80, 80);
			ctficon = new ImageIcon(ip.getImage());
		}
		return ctficon;
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
