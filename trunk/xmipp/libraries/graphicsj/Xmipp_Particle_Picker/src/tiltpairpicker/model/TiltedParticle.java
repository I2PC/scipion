package tiltpairpicker.model;

import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.MicrographParticle;

public class TiltedParticle extends MicrographParticle {
	
	private UntiltedParticle untiltedparticle;
	
	public TiltedParticle(int x, int y, UntiltedParticle untiltedparticle)
	{
		super(x, y, ((UntiltedMicrograph)untiltedparticle.getMicrograph()).getTiltedMicrograph());
		this.untiltedparticle = untiltedparticle;
	}
	
	public UntiltedParticle getUntiltedParticle()
	{
		return untiltedparticle;
	}

}
