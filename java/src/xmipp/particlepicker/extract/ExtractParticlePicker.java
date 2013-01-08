package xmipp.particlepicker.extract;

import java.util.ArrayList;
import java.util.List;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.particlepicker.IJCommand;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.ParticlePicker;
import xmipp.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.particlepicker.training.model.FamilyState;

public class ExtractParticlePicker extends ParticlePicker
{

	public static void main(String[] args)
	{
		ExtractParticlePicker.open(args[0]);
	}
	
	private ArrayList<ExtractMicrograph> micrographs;
	private ExtractMicrograph micrograph;

	public ExtractParticlePicker(String selfile, FamilyState mode)
	{
		super(selfile, mode);
		loadParticles();
		if(filters.isEmpty())//user just started manual mode and has no filter, I select gaussian blur by default, will be applied when window opens
			filters.add(new IJCommand("Gaussian Blur...", "sigma=2"));
	}

	@Override
	public void loadEmptyMicrographs()
	{
		micrographs = new ArrayList<ExtractMicrograph>();
		MetaData md = new MetaData(selfile);
		ExtractParticle p;
		String fileiter;
		boolean exists;
		ExtractMicrograph current = null;
		boolean existspsd = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
		boolean existsctf = md.containsLabel(MDLabel.MDL_CTF_MODEL);
		String psd = null, ctf = null;
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
				if (existspsd)
					psd = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_CTF_MODEL, id);
				current = new ExtractMicrograph(fileiter, psd, ctf);
				
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
			enabled = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1)? true: false;
			p = new ExtractParticle(x, y, current, enabled);
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
	public ExtractMicrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		micrograph = (ExtractMicrograph)m;

	}

	@Override
	public void saveConfig()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void loadConfig()
	{
		setMicrograph(micrographs.get(0));
	}

	public static void open(String filename)
	{
		ExtractParticlePicker picker = new ExtractParticlePicker(filename,  FamilyState.Extract);
		new ExtractPickerJFrame(picker);
	}

	public int getParticlesTotal()
	{
		int particles = 0;
		for(ExtractMicrograph m: micrographs)
			particles += m.getParticles().size();
		return particles;
	}

	public void resetAllMicrographs()
	{
		for (ExtractMicrograph m : micrographs)
			m.reset();
		setChanged(true);
		
	}
	
	
}
