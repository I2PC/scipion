package xmipp.viewer.particlepicker.tiltpair.model;

import xmipp.viewer.particlepicker.training.model.TrainingParticle;

public class TiltedParticle extends TrainingParticle {
	
	private UntiltedParticle untiltedparticle;
	
	public TiltedParticle(int x, int y, UntiltedParticle untiltedparticle)
	{
		super(x, y, untiltedparticle.getParticlePicker(), ((UntiltedMicrograph)untiltedparticle.getMicrograph()).getTiltedMicrograph());
		this.untiltedparticle = untiltedparticle;
	}
	
	public UntiltedParticle getUntiltedParticle()
	{
		return untiltedparticle;
	}



}
