package tiltpairpicker.model;

import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.Particle;

public class TiltedParticle extends Particle {
	
	private UntiltedParticle untiltedparticle;
	
	public TiltedParticle(int x, int y, TiltedMicrograph micrograph, UntiltedParticle untiltedparticle)
	{
		super(x, y, micrograph);
		this.untiltedparticle = untiltedparticle;
	}
	
	public UntiltedParticle getUntiltedParticle()
	{
		return untiltedparticle;
	}

}
