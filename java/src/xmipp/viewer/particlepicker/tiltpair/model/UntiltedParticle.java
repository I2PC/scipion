package xmipp.viewer.particlepicker.tiltpair.model;

import xmipp.viewer.particlepicker.ParticlePicker;
import xmipp.viewer.particlepicker.training.model.TrainingParticle;


public class UntiltedParticle extends TrainingParticle{
	
	private TiltedParticle tiltedparticle;
	private boolean added = false;
	
		 
	public boolean isAdded()
	{
		return added;
	}
	
	public UntiltedParticle(int x, int y, UntiltedMicrograph micrograph, ParticlePicker picker)
	{
		super(x, y, picker, micrograph);
		
	}
	
	public void setTiltedParticle(TiltedParticle tiltedparticle)
	{
		this.tiltedparticle = tiltedparticle;
	}
	
	public TiltedParticle getTiltedParticle()
	{
		return tiltedparticle;
	}

	public void setAdded(boolean added)
	{
		this.added = added;
	}
	
	
	

}
