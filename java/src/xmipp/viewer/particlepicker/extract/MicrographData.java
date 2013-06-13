package xmipp.viewer.particlepicker.extract;

import ij.gui.ImageWindow;

import java.util.ArrayList;
import java.util.List;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.Micrograph;

public class MicrographData
{
	private String mdfile;
	private List<Particle> particles;
	private Micrograph micrograph;
	private String micfile;

	public MicrographData(String mdfile, String micfile)
	{
		this.mdfile = mdfile;
		this.particles = new ArrayList<Particle>();
		this.micfile = micfile;
		this.micrograph = new Micrograph(micfile)
		{

			@Override
			public void reset()
			{
				// TODO Auto-generated method stub

			}

			@Override
			public boolean hasData()
			{
				return !particles.isEmpty();
			}
		};
		loadData();
	}

	public Micrograph getMicrograph()
	{
		return micrograph;
	}
	
	public static void openMicrographs(String mdfile)
	{
		
		
		List<MicrographData> micdatas = MicrographData.getMicrographData(mdfile);
		MicrographCanvas mc;
		for(MicrographData pl: micdatas)
		{
			mc = new MicrographCanvas(pl, 30);
			mc.display();
		}
	}
	


	private void loadData()
	{
		particles.clear();
		MetaData md = new MetaData(mdfile);
		Particle p;
		int x, y;
		String file;
		for (long id : md.findObjects())
		{
			file = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
			if (this.micfile.equals(file))
			{
				x = md.getValueInt(MDLabel.MDL_XCOOR, id);
				y = md.getValueInt(MDLabel.MDL_YCOOR, id);
				p = new Particle(x, y);
				particles.add(p);
			}
		}
	}

	public List<Particle> getParticles()
	{
		return particles;
	}

	public String getMicrographFile()
	{
		return micfile;
	}

	public static List<MicrographData> getMicrographData(String file)
	{
		List<MicrographData> loaders = new ArrayList<MicrographData>();
		MetaData md = new MetaData(file);
		Particle p;
		int x, y;
		String micfileiter;
		boolean exists;
		MicrographData data = null;
		for (long id : md.findObjects())
		{
			exists = false;
			micfileiter = md.getValueString(MDLabel.MDL_MICROGRAPH, id);
			for (MicrographData dataiter : loaders)
				if (dataiter.getMicrographFile().equals(micfileiter))
				{
					exists = true;
					data = dataiter;
					break;
				}
			if (!exists)
			{
				data = new MicrographData(file, micfileiter);
				loaders.add(data);
			}
			x = md.getValueInt(MDLabel.MDL_XCOOR, id);
			y = md.getValueInt(MDLabel.MDL_YCOOR, id);
			p = new Particle(x, y);
			data.addParticle(p);
		}
		return loaders;
	}

	private void addParticle(Particle p)
	{
		particles.add(p);
	}

}
