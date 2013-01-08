package xmipp.particlepicker.particles;

import java.util.ArrayList;
import java.util.List;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePicker;
import xmipp.particlepicker.training.model.FamilyState;

public class ExtractParticlePicker extends ParticlePicker
{

	private ArrayList<ExtractMicrograph> micrographs;
	private Micrograph micrograph;

	public ExtractParticlePicker(String selfile, String outputdir, FamilyState mode)
	{
		super(selfile, outputdir, mode);
		loadParticles();
	}

	@Override
	public void loadEmptyMicrographs()
	{
		micrographs = new ArrayList<ExtractMicrograph>();
		MetaData md = new MetaData(selfile);
		ExtractParticle p;
		int x, y;
		String fileiter;
		boolean exists, enabled;
		ExtractMicrograph current = null;
		for (long id : md.findObjects())
		{
			exists = false;
			fileiter = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
			for (ExtractMicrograph iter : micrographs)
				if (iter.getFile().equals(fileiter))
				{
					exists = true;
					current = iter;
					break;
				}
			if (!exists)
			{
				current = new ExtractMicrograph(fileiter);
				micrographs.add(current);
			}

		}

	}
	
	public void loadParticles()
	{
		MetaData md = new MetaData(selfile);
		ExtractParticle p;
		int x, y;
		String fileiter;
		boolean enabled;
		ExtractMicrograph current = null;
		for (long id : md.findObjects())
		{
			fileiter = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
			for (ExtractMicrograph iter : micrographs)
				if (iter.getFile().equals(fileiter))
				{
			
					current = iter;
					break;
				}
			
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			enabled = md.getValueBoolean(MDLabel.MDL_ENABLED, id);
			p = new ExtractParticle(x, y, enabled);
			current.addParticle(p);

		}
	}

	@Override
	public List<? extends Micrograph> getMicrographs()
	{
		return micrographs;
	}

	@Override
	public void saveData(Micrograph m)
	{
	
	}

	



	@Override
	public Micrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		micrograph = m;

	}

	@Override
	public void saveConfig()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void loadConfig()
	{
		// TODO Auto-generated method stub
		
	}

	public static void open(String filename)
	{
		// TODO Auto-generated method stub
		
	}

}
