package ppicker;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import xmipp.MDLabel;
import xmipp.MetaData;

public class PPData {
	
	private List<Family> families;
	
	public PPData()
	{
		this.families = new ArrayList<Family>();
		families.add(Family.getDefaultFamily());
	}
	
	
	public void saveData(String filename)
	{
		long id;
		try {
			MetaData md = new MetaData();
			for(Family f: families)
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, f.getName(), id);
				md.setValueInt(MDLabel.MDL_XINT, f.getColor().getRGB(), id);
				md.setValueInt(MDLabel.MDL_YINT, f.getRadius(), id);
			}
			md.writeBlock("groups@" + filename);
			for(Family g: families)
			{
				md = new MetaData();
				for (Particle p: g.getParticles()) {
					id = md.addObject();
					md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
					md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
				}
				md.writeBlock(g.getName() + "@" + filename);
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public List<Family> getFamilies() {
		return families;
	}

	

	public void loadData(String filename)
	{
		families.clear();
		Family family;
		int rgb, radius;
		String gname;		
		try {
			MetaData md = new MetaData("groups@" + filename);
			long[] ids = md.findObjects();
			for (long id: ids) {				
				gname = md.getValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, id);
				rgb = md.getValueInt(MDLabel.MDL_XINT, id);
				radius = md.getValueInt(MDLabel.MDL_YINT, id);
				family = new Family(gname, new Color(rgb), radius);
				fillFamilyParticles(filename, family);
				families.add(family);
			}				
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void fillFamilyParticles(String filename, Family family)
	{
		int x, y;
		MetaData md;
		try {
			md = new MetaData(family.getName() + "@" + filename);
			long[] ids = md.findObjects();
			for (long id: ids) {				
				x = md.getValueInt(MDLabel.MDL_XINT, id);
				y = md.getValueInt(MDLabel.MDL_YINT, id);
				family.addParticle(new Particle(x, y, family));
			}		
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public boolean existsGroupName(String name) {
		for (Family g2 : getFamilies())
			if (g2.getName().equalsIgnoreCase(name))
				return true;
		return false;
	}

}
