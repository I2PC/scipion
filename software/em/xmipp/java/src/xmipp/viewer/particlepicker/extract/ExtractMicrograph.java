package xmipp.viewer.particlepicker.extract;

import java.util.ArrayList;
import java.util.List;

import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.CtfInfo;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.PickerParticle;

public class ExtractMicrograph extends Micrograph
{

	private ArrayList<ExtractParticle> particles;

	public ExtractMicrograph(String file, CtfInfo ctfInfo)
	{
		super(file, ctfInfo);
		particles = new ArrayList<ExtractParticle>();
	}

	@Override
	public boolean hasData()
	{
		return !particles.isEmpty();
	}

    @Override
    public List<? extends PickerParticle> getParticleList() {
        return getParticles();
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
