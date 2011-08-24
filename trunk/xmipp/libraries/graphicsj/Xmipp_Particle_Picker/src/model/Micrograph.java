package model;

import ij.ImagePlus;

import java.awt.Image;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;

public class Micrograph {
	
	private String filename;
	private String name;
	private ImagePlus image;
	private String outputfname;
	private static String ext = ".pos";
	private String ctf;
	private ImageIcon ctficon;
	private String aoutputfname;
	private boolean autopicking = false;
	private List<MicrographFamilyData> mfdatas;
	
	
	public Micrograph(String filename, String ctf) {
		this.filename = filename;
		this.name = getName(filename);
		mfdatas = new ArrayList<MicrographFamilyData>();
		this.outputfname = ExecutionEnvironment.getOutputPath(name + ext);
		this.aoutputfname = ExecutionEnvironment.getOutputPath(name + "_auto" + ext);
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
		return ExecutionEnvironment.getMicrographsSelFile();
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
			Image i = ip.getImage().getScaledInstance(150, 150, Image.SCALE_SMOOTH);
			ctficon = new ImageIcon(i);
			
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
	
	public List<MicrographFamilyData> getFamiliesData()
	{
		return mfdatas;
	}

	public Particle getParticle(int x, int y)
	{
		for(MicrographFamilyData mfd: mfdatas)
			for(Particle p: mfd.getParticles())
			{
				if (p.contains(x, y)) 
				return p;
			}
			return null;
	}
	
	public MicrographFamilyData getFamilyData(Family f)
	{
		for(MicrographFamilyData fp: mfdatas)
			if(f.equals(fp.getFamily()))
				return fp;
		MicrographFamilyData mfd = new MicrographFamilyData(f); 
		mfdatas.add(mfd);
		return mfd;
	}

	public void addParticle(Particle p)
	{
		getFamilyData(p.getFamily()).addParticle(p);
	}
	
	
	public void removeParticle(Particle p)
	{
		getFamilyData(p.getFamily()).removeParticle(p);
	}
	
	public String toString()
	{
		return name;
	}

	

}
