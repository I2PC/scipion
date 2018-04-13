package xmipp.viewer.particlepicker.extract;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import xmipp.ij.commons.IJCommand;
import xmipp.ij.commons.XmippImageJ;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.viewer.particlepicker.*;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.windows.GalleryJFrame;

public class ExtractParticlePicker extends ParticlePicker
{

	private ArrayList<ExtractMicrograph> micrographs;
	private ExtractMicrograph micrograph;
	private ArrayList<ColorHelper> colorby;


	public ExtractParticlePicker(String selfile, int size, ParticlePickerParams params)
	{
		super(selfile, params);
		setSize(size);
		loadParticles();
		if (filters.isEmpty())
			filters.add(new IJCommand(XmippImageJ.gaussianBlurFilter, "sigma=2"));
	}
	
	public ExtractParticlePicker(String block, String selfile, int size, ParticlePickerParams params)
	{
		super(block, selfile, ".", params);
		setSize(size);
		loadParticles();
		if (filters.isEmpty())
			filters.add(new IJCommand(XmippImageJ.gaussianBlurFilter, "sigma=2"));
	}
	
	

	@Override
	public void loadEmptyMicrographs()
	{
		micrographs = new ArrayList<ExtractMicrograph>();
		String path = (block == null)? selfile: block + "@" + selfile;
		MetaData md = new MetaData(path);
		String fileiter;
		boolean exists;
		ExtractMicrograph current = null;
		boolean existspsd = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
		boolean existsctf = md.containsLabel(MDLabel.MDL_CTF_MODEL);
		String psd = null, ctf = null;
        CtfInfo ctfInfo;

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
			    ctfInfo = (existsctf || existspsd) ? new CtfInfo() : null;
				if (existspsd)
					ctfInfo.psd = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				if (existsctf)
					ctfInfo.ctf = md.getValueString(MDLabel.MDL_CTF_MODEL, id);
				current = new ExtractMicrograph(fileiter, ctfInfo);

				micrographs.add(current);
			}

		}
		loadColumns(md);
		md.destroy();

	}
	
	public void loadColumns(MetaData md)
	{
		colorby = new ArrayList<ColorHelper>();
                
		loadColumn(MDLabel.MDL_ZSCORE, "ZScore", md);
		loadColumn(MDLabel.MDL_ZSCORE_SHAPE1, "ZScore-Shape1", md);
		loadColumn(MDLabel.MDL_ZSCORE_SHAPE2, "ZScore-Shape2", md);
		loadColumn(MDLabel.MDL_ZSCORE_SNR1, "ZScore-SNR1", md);
		loadColumn(MDLabel.MDL_ZSCORE_SNR2, "ZScore-SNR2", md);
		loadColumn(MDLabel.MDL_ZSCORE_HISTOGRAM, "ZScore-Hist", md);
                if(colorby.isEmpty())
                    throw new IllegalArgumentException("No score available to sort particles");
		
	}
	
	public void loadColumn(int column, String name, MetaData md)
	{
		boolean exists = md.containsLabel(column);
		if(exists)
		{
			colorby.add(new ColorHelper(column, name, md));
		}
	}
	
	

	public void loadParticles()
	{
		String path = (block == null)? selfile: block + "@" + selfile;
		MetaData md = new MetaData(path);
		ExtractParticle p;
		int x, y;
		String fileiter;
		boolean enabled;
		ExtractMicrograph current = null;
		double zscore, zscore_shape1, zscore_shape2, zscore_snr1, zscore_snr2, zscore_hist;
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
			zscore_shape1 = md.getValueDouble(MDLabel.MDL_ZSCORE_SHAPE1, id);
			zscore_shape2 = md.getValueDouble(MDLabel.MDL_ZSCORE_SHAPE2, id);
			zscore_snr1 = md.getValueDouble(MDLabel.MDL_ZSCORE_SNR1, id);
			zscore_snr2 = md.getValueDouble(MDLabel.MDL_ZSCORE_SNR2, id);
			zscore_hist = md.getValueDouble(MDLabel.MDL_ZSCORE_HISTOGRAM, id);
			p = new ExtractParticle(id, x, y, current, enabled, zscore, zscore_shape1, zscore_shape2, zscore_snr1, zscore_snr2, zscore_hist);
			current.addParticle(p);

		}
		md.destroy();
	}

	@Override
	public List<ExtractMicrograph> getMicrographs()
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
			MetaData md = new MetaData();
			
			for (ExtractParticle p : micrograph.getParticles())
			{
				id = p.getId();
				md.setValueInt(MDLabel.MDL_XCOOR, p.getX(), id);
				md.setValueInt(MDLabel.MDL_YCOOR, p.getY(), id);
				md.setValueInt(MDLabel.MDL_ENABLED, p.isEnabled()? 1: -1, id);
			}
			String path = (block == null)? selfile: block + "@" + selfile;
			md.write(path);
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

	public static ExtractPickerJFrame open(String block, String filename, int size, GalleryJFrame galleryfr)
	{
		try
		
		{
			ExtractParticlePicker picker = new ExtractParticlePicker(block, filename, size, new ParticlePickerParams(Mode.Extract));
			return new ExtractPickerJFrame(picker, galleryfr);
		}
		catch(Exception e)
		{
			XmippDialog.showError(null,e.getMessage());
			return null;
		}
	}

        @Override
	public int getParticlesCount()
	{
		int particles = 0;
		for (ExtractMicrograph m : micrographs)
			particles += m.getParticles().size();
		return particles;
	}

	

	public ColorHelper[] getColumns()
	{
		return colorby.toArray(new ColorHelper[]{});
	}

	@Override
	public boolean isValidSize(JFrame parent, int size)
	{
		for (ExtractParticle p : getMicrograph().getParticles())
			if (!getMicrograph().fits(p.getX(), p.getY(), size))
				return false;
		return true;
	}

	

  

}
