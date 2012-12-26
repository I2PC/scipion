package xmipp.viewer.windows;

import java.util.ArrayList;
import java.util.List;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.jni.Particle;
import xmipp.particlepicker.Micrograph;

public class ParticlesLoader
{
	private String mdfile;
	private List<Particle> particles;
	private Micrograph micrograph;
	private String micfile;

	public ParticlesLoader(String mdfile, String micfile)
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

}
