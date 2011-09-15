package tiltpairpicker.model;

import ij.ImagePlus;
import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.Particle;


public class UntiltedParticle extends Particle{
	
	private TiltedParticle tiltedparticle;
	
	
	
	public UntiltedParticle(int x, int y, UntiltedMicrograph micrograph)
	{
		super(x, y, micrograph);
	}
	
	public void setTiltedParticle(TiltedParticle tiltedparticle)
	{
		this.tiltedparticle = tiltedparticle;
	}
	

	

}
