package xmipp.viewer.particlepicker.training.model;


import ij.ImagePlus;

import java.awt.Image;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import xmipp.ij.commons.XmippIJUtil;
import xmipp.jni.Filename;
import xmipp.viewer.particlepicker.Family;
import xmipp.viewer.particlepicker.Micrograph;

public class TrainingMicrograph extends Micrograph{
	
	private String psd;
	private ImageIcon ctficon;
	private boolean autopicking = false;
	private List<MicrographFamilyData> mfdatas;
	private String autofilename;
	private String ctf;
	
	public TrainingMicrograph(String filename, String psd, String ctf, List<Family> families, FamilyState mode) {
		this(filename, psd, ctf, families, new ArrayList<MicrographFamilyData>(), mode);
	}
	

	public TrainingMicrograph(String file, String psd, String ctf, List<Family> families, List<MicrographFamilyData> mfd, FamilyState mode) {
		super(file);
		this.psd = psd;
		this.ctf = ctf;
		mfdatas = mfd;
		autofilename = getName() + "_auto" + ext;
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
	

	
	public String getAutoPosFile()
	{
		return autofilename;
	}
	
	public Icon getCTFIcon()
	{
		String file;
		if(ctficon == null)
		{
			if(psd == null || !(new File(psd).exists()))
				file = (Filename.getXmippPath("resources" + File.separator + "no-image.jpg"));
			else
				file = psd;
			Image image = XmippIJUtil.getImagePlus(file).getImage().getScaledInstance(120, 110, Image.SCALE_SMOOTH);
			ctficon = new ImageIcon(image);
			
		}
		return ctficon;
	}
	

	public ImagePlus getPSDImage()
	{
			if(psd == null || !(new File(psd).exists()))
				return null;
			return XmippIJUtil.getImagePlus(psd);
			
	}
	
	public String getPSD()
	{
		return psd;
	}
	
	public String getCTF()
	{
		return ctf;
	}
	
	
	public List<MicrographFamilyData> getFamiliesData()
	{
		return mfdatas;
	}

	public TrainingParticle getParticle(int x, int y)
	{
		for(MicrographFamilyData mfd: mfdatas)
		{
			for(TrainingParticle p: mfd.getManualParticles())
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
		MicrographFamilyData mfd = new MicrographFamilyData(this, f); 		
		mfdatas.add(mfd);	
		return mfd;
	}
	
	
	public void addManualParticle(TrainingParticle p, TrainingPicker ppicker, boolean center, boolean totemplates)
	{
		getFamilyData(p.getFamily()).addManualParticle(p, ppicker, center, totemplates);
		
	}
	
	
	public void removeParticle(TrainingParticle p, TrainingPicker ppicker)
	{
		TrainingParticle tp = (TrainingParticle)p;
		MicrographFamilyData mfd = getFamilyData(ppicker.getFamily());
		mfd.removeParticle(tp, ppicker);
		
	}
	
	public void addAutomaticParticle(AutomaticParticle p)
	{
		addAutomaticParticle(p, false);
	}
	
	public void addAutomaticParticle(AutomaticParticle p, boolean imported)
	{
		getFamilyData(p.getFamily()).addAutomaticParticle(p, imported);
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
//			if(mfd.getState() != MicrographFamilyState.Available)
			if(!mfd.isEmpty())
				return true;
		return false;
	}
	
	public void reset()
	{
		for(MicrographFamilyData mfd: mfdatas)
			if(mfd != null)
				mfd.reset();
	}


	public void removeFamilyData(Family family) {
		mfdatas.remove(getFamilyData(family));
		
	}

	public void removeParticles(int x, int y, TrainingPicker ppicker)
	{
		List<TrainingParticle> particles = new ArrayList<TrainingParticle>();
		for(MicrographFamilyData mfd: mfdatas)
		{
			for(TrainingParticle p: mfd.getManualParticles())
				if (p.contains(x, y)) 
					particles.add(p);
			for(TrainingParticle p: particles)
				removeParticle(p, ppicker);
			particles.clear();
			for(AutomaticParticle p: mfd.getAutomaticParticles())
				if (p.contains(x, y)) 
					particles.add(p);

			for(TrainingParticle p: particles)
				removeParticle(p, ppicker);
		}
		
	}
	
	


	
	


}
