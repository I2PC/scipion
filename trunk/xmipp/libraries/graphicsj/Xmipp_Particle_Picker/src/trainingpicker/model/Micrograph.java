package trainingpicker.model;

import ij.ImagePlus;

import java.awt.Image;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;

public class Micrograph {
	
	private String file;
	private String name;
	private ImagePlus image;
	private String outputfilename;
	private static String ext = ".pos";
	private String ctf;
	private ImageIcon ctficon;
	private boolean autopicking = false;
	private List<MicrographFamilyData> mfdatas;
	private String autofilename;
	
	public Micrograph(String filename, String ctf, List<Family> families, FamilyState mode) {
		this(filename, ctf, families, new ArrayList<MicrographFamilyData>(), mode);
	}
	

	public Micrograph(String file, String ctf, List<Family> families, List<MicrographFamilyData> mfd, FamilyState mode) {
		
		this.file = file;
		if(!new File(file).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("file", file));
		this.name = getName(file);
		this.ctf = ctf;
		mfdatas = mfd;
		this.outputfilename = name + ext;
		autofilename = name + "_auto" + ext;
		MicrographFamilyState state = (mode == FamilyState.Review)? MicrographFamilyState.Review : MicrographFamilyState.Available;
		for(Family f: families)
			mfdatas.add(new MicrographFamilyData(this, f, state));

	}
	

	
	
	void setFamiliesState(List<MicrographFamilyData> mfdatas)
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
	
	public static String getName(String file)
	{
		String[] tokens = file.split(File.separator);
		if(tokens.length < 2)
			throw new IllegalArgumentException("Name for micrograph" +
					"is taken from parent dir, invalid path " + file);
		return  tokens[tokens.length - 2];
	}
	

	
	public String getAutoFilename()
	{
		return autofilename;
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
	
	public Icon getCTFIcon()
	{
		String file;
		if(ctficon == null)
		{
			if(ctf == null || !(new File(ctf).exists()))
				file = (ParticlePicker.getXmippPath("resources" + File.separator + "no-image.jpg"));
			else
				file = ctf;
			Image image = new ImagePlus(file).getImage().getScaledInstance(120, 110, Image.SCALE_SMOOTH);
			ctficon = new ImageIcon(image);
			
		}
		return ctficon;
	}
	
	
	public String getFile() {
		return file;
	}
	public void setFile(String file) {
		this.file = file;
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
			
		}
		return null;
	}
	
	public AutomaticParticle getAutomaticParticle(int x, int y, double threshold)
	{
		for(MicrographFamilyData mfd: mfdatas)
			for(AutomaticParticle p: mfd.getAutomaticParticles())
				if (!p.isDeleted() && p.getCost() >= threshold && p.contains(x, y)) 
					return p;
		return null;
	}
	
	public MicrographFamilyData getFamilyData(Family f)
	{
		for(MicrographFamilyData fp: mfdatas)
			if(f.equals(fp.getFamily()))
				return fp;
		return null;
	}

	public void addManualParticle(Particle p)
	{
		getFamilyData(p.getFamily()).addManualParticle(p);
	}
	
	
	public void removeParticle(Particle p, ParticlePicker ppicker)
	{
		getFamilyData(p.getFamily()).removeParticle(p, ppicker);
	}
	
	public void addAutomaticParticle(AutomaticParticle p)
	{
		getFamilyData(p.getFamily()).addAutomaticParticle(p);
	}
	
	
	public String toString()
	{
		return name;
	}

	public boolean hasManualParticles()
	{
		for(MicrographFamilyData fp: mfdatas)
			if(fp != null && fp.getManualParticles().size() > 0)
				return true;
		return false;
	}
	
	public boolean hasAutomaticParticles()
	{
		for(MicrographFamilyData fp: mfdatas)
			if(fp != null && !fp.getAutomaticParticles().isEmpty())
				return true;
		return false;
	}
	
	public boolean isPickingAvailable(Family f) {
			return getFamilyData(f).isPickingAvailable();
	}
	
	public boolean hasData()
	{
		for(MicrographFamilyData mfd: mfdatas)
			if(mfd.getState() != MicrographFamilyState.Available)
				return true;
		return false;
	}
	
	public void reset()
	{
		for(MicrographFamilyData mfd: mfdatas)
			if(mfd != null)
				mfd.reset();
	}


}
