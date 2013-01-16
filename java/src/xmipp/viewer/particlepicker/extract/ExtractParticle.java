package xmipp.viewer.particlepicker.extract;

import xmipp.jni.Particle;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.PickerParticle;

public class ExtractParticle extends PickerParticle
{
	private boolean enabled;
	private long id;

	public ExtractParticle(long id, int x, int y, Micrograph m, boolean enabled)
	{
		super(x, y, m);
		this.enabled = enabled;
		this.id = id;
	}
	
	public long getId()
	{
		return id;
	}

	public ExtractParticle(long id, int x, int y, Micrograph m)
	{
		this(id, x, y, m, true);
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
