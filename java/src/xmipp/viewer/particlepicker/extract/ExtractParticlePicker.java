package xmipp.viewer.particlepicker.extract;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.viewer.particlepicker.IJCommand;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;
import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyData;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;

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
		if (filters.isEmpty())// user just started manual mode and has no
								// filter, I select gaussian blur by default,
								// will be applied when window opens
			filters.add(new IJCommand("Gaussian Blur...", "sigma=2"));
	}

	@Override
	public void loadEmptyMicrographs()
	{
		micrographs = new ArrayList<ExtractMicrograph>();
		MetaData md = new MetaData(selfile);
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
		double zscore, zscore_shape, zscore_snr1, zscore_snr2, zscore_hist;
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
			enabled = (md.getValueInt(MDLabel.MDL_ENABLED, id) == 1) ? true : false;
			zscore = md.getValueDouble(MDLabel.MDL_ZSCORE, id);
			zscore_shape = md.getValueDouble(MDLabel.MDL_ZSCORE_SHAPE, id);
			zscore_snr1 = md.getValueDouble(MDLabel.MDL_ZSCORE_SNR1, id);
			zscore_snr2 = md.getValueDouble(MDLabel.MDL_ZSCORE_SNR2, id);
			zscore_hist = md.getValueDouble(MDLabel.MDL_ZSCORE_HISTOGRAM, id);
			p = new ExtractParticle(id, x, y, current, enabled, zscore, zscore_shape, zscore_snr1, zscore_snr2, zscore_hist);
			current.addParticle(p);

		}
	}

	@Override
	public List<? extends Micrograph> getMicrographs()
	{
		return micrographs;
	}

	public void saveData()
	{
		if (isChanged())
			for(ExtractMicrograph m: micrographs)
				saveData(m);
	}
	
	@Override
	public void saveData(Micrograph m)
	{
		long id;
		micrograph = (ExtractMicrograph)m;
		try
		{
			MetaData md = new MetaData(selfile);
			
			for (ExtractParticle p : micrograph.getParticles())
			{
				id = p.getId();
				md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
				md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
				md.setValueInt(MDLabel.MDL_ENABLED, p.isEnabled()? 1: -1, id);
			}
			md.write(selfile);
			md.destroy();

		}
		catch (Exception e)
		{
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	@Override
	public ExtractMicrograph getMicrograph()
	{
		return micrograph;
	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		micrograph = (ExtractMicrograph) m;

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
		ExtractParticlePicker picker = new ExtractParticlePicker(filename, FamilyState.Extract);
		new ExtractPickerJFrame(picker);
	}

	public int getParticlesTotal()
	{
		int particles = 0;
		for (ExtractMicrograph m : micrographs)
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
