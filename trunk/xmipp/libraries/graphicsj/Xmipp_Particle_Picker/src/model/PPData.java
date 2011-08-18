package model;

import gui.PPCanvas;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import xmipp.MDLabel;
import xmipp.MetaData;

public class PPData {
	
	private List<Family> families;
	private List<Micrograph> micrographs;
	private static PPData ppdata;
	
	private PPData()
	{
		this.families = new ArrayList<Family>();
		this.micrographs = new ArrayList<Micrograph>();
		loadFamilyData();
		loadMicrographsData();
	}
	
	public static  PPData getInstance()
	{
		if(ppdata == null)
			ppdata = new PPData();
		return ppdata;
	}
	
	public List<Family> getFamilies() {
		return families;
	}
	
	public List<Micrograph> getMicrographs()
	{
		return micrographs;
	}

	
	public void saveFamilyData()
	{
		long id;
		String filename = PPConfiguration.getFamiliesXMD();
		try {
			MetaData md = new MetaData();
			for(Family f: families)
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, f.getName(), id);
				md.setValueInt(MDLabel.MDL_XINT, f.getColor().getRGB(), id);
				md.setValueInt(MDLabel.MDL_YINT, f.getSize(), id);
			}
			md.write("families@" + filename);
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
		
	}
	


	public void loadFamilyData()
	{
		families.clear();
		String xmd = PPConfiguration.getFamiliesXMD();
		if(!new File(xmd).exists())
		{
			families.add(Family.getDefaultFamily());
			return;
		}
		
		Family family;
		int rgb, size;
		String gname;		
		try {
			MetaData md = new MetaData("families@" + xmd);
			long[] ids = md.findObjects();
			for (long id: ids) {				
				gname = md.getValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, id);
				rgb = md.getValueInt(MDLabel.MDL_XINT, id);
				size = md.getValueInt(MDLabel.MDL_YINT, id);
				family = new Family(gname, new Color(rgb), size);
				families.add(family);
			}				
			if(families.size() == 0)
				throw new IllegalArgumentException(String.format("No families specified on %s", xmd));
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	
	public Family getFamily(String name)
	{
		for (Family f : getFamilies())
			if (f.getName().equalsIgnoreCase(name))
				return f;
		return null;
	}
	
	public boolean existsFamilyName(String name) {
		return getFamily(name)!= null;
	}
	
	
	public void loadMicrographsData()
	{
		String xmd = PPConfiguration.getMicrographsXMD();
		micrographs.clear();
		Micrograph micrograph;
		String name, filename;		
		try {
			MetaData md = new MetaData(xmd);
			int count = 1; 
			long[] ids = md.findObjects();
			for (long id: ids) {
				
				filename = PPConfiguration.getMicrographPath(md.getValueString(MDLabel.MDL_IMAGE, id));
				name = "micrograph" + count;
				micrograph = new Micrograph(filename, name);
				loadParticles(micrograph);
				micrographs.add(micrograph);
				count ++;
			}
			if(micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", xmd));
			
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
		
	}
	
	public void saveParticles(Micrograph micrograph) {
		long id;
		try {

				MetaData md = new MetaData();
				for (Particle p: micrograph.getParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					md.setValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, p.getFamily().getName(), id);
				}
				md.write("particles@" + micrograph.getXMD());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}

	public void loadParticles(Micrograph micrograph) {
		try {
			int x, y;
			String fname;
			Family family;
			Particle particle;
			if(!new File(micrograph.getXMD()).exists())
				return;
			MetaData md = new MetaData("particles@" + micrograph.getXMD());
			long[] ids = md.findObjects();
			for (long id: ids) {				
				
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				fname = md.getValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, id);
				family = getFamily(fname);
				particle = new Particle(x, y, family, micrograph);
				micrograph.addParticle(particle);
			}				
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
		
	}

}
