package xmipp.particlepicker.particles;

import xmipp.jni.Particle;

public class ExtractParticle extends Particle
{
	private boolean enabled;

	public ExtractParticle(int x, int y, boolean enabled)
	{
		super(x, y);
		this.enabled = enabled;
	}
	
	public ExtractParticle(int x, int y)
	{
		super(x, y);
		this.enabled = true;
	}

	public boolean isEnabled()
	{
		return enabled;
	}

	public void setEnabled(boolean enabled)
	{
		this.enabled = enabled;
	}
	
	

}
