package tiltpairpicker.model;

import trainingpicker.model.Family;
import trainingpicker.model.TrainingParticle;


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
