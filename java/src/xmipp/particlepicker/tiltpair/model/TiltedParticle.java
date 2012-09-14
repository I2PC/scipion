package xmipp.particlepicker.tiltpair.model;

import xmipp.particlepicker.training.model.TrainingParticle;

public class TiltedParticle extends TrainingParticle {
	
	private UntiltedParticle untiltedparticle;
	
	public TiltedParticle(int x, int y, UntiltedParticle untiltedparticle)
	{
		super(x, y, untiltedparticle.getFamily(), ((UntiltedMicrograph)untiltedparticle.getMicrograph()).getTiltedMicrograph());
		this.untiltedparticle = untiltedparticle;
	}
	
	public UntiltedParticle getUntiltedParticle()
	{
		return untiltedparticle;
	}

}
