package xmipp.viewer.particlepicker.extract;

import java.util.ArrayList;
import java.util.List;

import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;

public class ExtractMicrograph extends Micrograph
{

	private ArrayList<ExtractParticle> particles;

	public ExtractMicrograph(String file, String psd, String ctf)
	{
		super(file, psd, ctf);
		particles = new ArrayList<ExtractParticle>();
	}

	@Override
	public boolean hasData()
	{
		return !particles.isEmpty();
	}

	@Override
	public void reset()
	{
		particles.clear();

	}

	public void addParticle(ExtractParticle p)
	{
		particles.add(p);

	}

	public List<ExtractParticle> getParticles()
	{
		return particles;
	}

	public void removeParticles(int x, int y, ParticlePicker picker)
	{
		List<ExtractParticle> removelist = new ArrayList<ExtractParticle>();

		for (ExtractParticle p : particles)
			if (p.contains(x, y, picker.getSize()))
				removelist.add(p);
		for (ExtractParticle p : removelist)
			removeParticle(p);
		

	}

	public void removeParticle(ExtractParticle particle)
	{
		particle.setEnabled(false);
	}

	public ExtractParticle getParticle(int x, int y, int size)
	{
		for(ExtractParticle p: particles)
			if(p.contains(x, y, size))
				return p;
		return null;
	}

}
