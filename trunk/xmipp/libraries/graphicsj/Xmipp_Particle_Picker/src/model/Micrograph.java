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
		this.outputfname = ParticlePicker.getOutputPath(name + ext);
		this.aoutputfname = ParticlePicker.getOutputPath(name + "_auto" + ext);
		this.ctf = ctf;
	}
	

	public Micrograph(String filename, String ctf, List<MicrographFamilyData> mfd) {
		this.filename = filename;
		this.name = getName(filename);
		mfdatas = mfd;
		this.outputfname = ParticlePicker.getOutputPath(name + ext);
		this.aoutputfname = ParticlePicker.getOutputPath(name + "_auto" + ext);
		this.ctf = ctf;
	}
	
	void setFamiliesData(List<MicrographFamilyData> mfdatas)
	{
		this.mfdatas = mfdatas;
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
		return ParticlePicker.getMicrographsSelFile();
	}
	
	public String getOFilename()
	{
		return outputfname;
	}
	
	public String getAutoOFilename()
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
		String file;
		if(ctficon == null)
		{
			if(ctf == null || !(new File(ctf).exists()))
				file = (ParticlePicker.getXmippPath("resources" + File.separator + "no-image.jpg"));
			else
				file = ctf;
			Image image = new ImagePlus(file).getImage().getScaledInstance(110, 110, Image.SCALE_SMOOTH);
			ctficon = new ImageIcon(image);
			
		}
		return ctficon;
	}
	
	public String getOutputRoot()
	{
		return ParticlePicker.getOutputPath(name);
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
		{
			for(Particle p: mfd.getManualParticles())
				if (p.contains(x, y)) 
					return p;
			for(AutomaticParticle p: mfd.getAutomaticParticles())
				if (p.contains(x, y) && !p.isDeleted() && p.getCost() >= 0) 
					return p;
		}
		return null;
	}
	
	public MicrographFamilyData getFamilyData(Family f)
	{
		for(MicrographFamilyData fp: mfdatas)
			if(f.equals(fp.getFamily()))
				return fp;
		MicrographFamilyData mfd = new MicrographFamilyData(this, f); 
		mfdatas.add(mfd);
		return mfd;
	}

	public void addManualParticle(Particle p)
	{
		getFamilyData(p.getFamily()).addManualParticle(p);
	}
	
	
	public void removeParticle(Particle p)
	{
		getFamilyData(p.getFamily()).removeParticle(p);
	}
	
	public void addAutomaticParticle(AutomaticParticle p)
	{
		getFamilyData(p.getFamily()).addAutomaticParticle(p);
	}
	
	
	public String toString()
	{
		return name;
	}

	public boolean isEmpty()
	{
		for(MicrographFamilyData fp: mfdatas)
			if(fp != null && fp.getManualParticles().size() > 0)
				return false;
		return true;
	}
	
	public boolean isPickingAvailable(Family f) {
			return getFamilyData(f).isPickingAvailable();
	}
	
	


}
