package xmipp.particlepicker.extract;

import xmipp.jni.Particle;
import xmipp.particlepicker.Micrograph;
import xmipp.particlepicker.PickerParticle;

public class ExtractParticle extends PickerParticle
{
	private boolean enabled;

	public ExtractParticle(int x, int y, Micrograph m, boolean enabled)
	{
		super(x, y, m);
		this.enabled = enabled;
	}
	
	public ExtractParticle(int x, int y, Micrograph m)
	{
		this(x, y, m, true);
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
