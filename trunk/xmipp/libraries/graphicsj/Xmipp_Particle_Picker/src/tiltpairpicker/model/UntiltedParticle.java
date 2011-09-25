package tiltpairpicker.model;

import ij.ImagePlus;
import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.MicrographParticle;


public class UntiltedParticle extends MicrographParticle{
	
	private TiltedParticle tiltedparticle;
	private boolean added = false;
	

	 
	public boolean isAdded()
	{
		return added;
	}
	
	public UntiltedParticle(int x, int y, UntiltedMicrograph micrograph)
	{
		super(x, y, micrograph);
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
