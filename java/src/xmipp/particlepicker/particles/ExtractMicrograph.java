package xmipp.particlepicker.particles;

import java.util.ArrayList;
import java.util.List;
import xmipp.particlepicker.Micrograph;

public class ExtractMicrograph extends Micrograph
{

	private ArrayList<ExtractParticle> particles;

	public ExtractMicrograph(String file)
	{
		super(file);
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

}
