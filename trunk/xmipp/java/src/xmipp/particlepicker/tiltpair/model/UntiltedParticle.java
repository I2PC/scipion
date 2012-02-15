package xmipp.particlepicker.tiltpair.model;

import xmipp.jni.Particle;
import xmipp.particlepicker.Family;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.XmippMessage;


public class UntiltedParticle extends TrainingParticle{
	
	private TiltedParticle tiltedparticle;
	private boolean added = false;
	
		 
	public boolean isAdded()
	{
		return added;
	}
	
	public UntiltedParticle(int x, int y, UntiltedMicrograph micrograph, Family family)
	{
		super(x, y, family, micrograph);
		if(!Particle.boxContainedOnImage(x, y, family.getSize(), micrograph.getTiltedMicrograph().getImagePlus()))
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("Tilted Pair Coordinates"));
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
