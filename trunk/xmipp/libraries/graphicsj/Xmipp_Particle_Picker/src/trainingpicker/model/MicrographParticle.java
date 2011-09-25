package trainingpicker.model;

import xmipp.Particle;



public class MicrographParticle extends Particle  {
	protected Micrograph micrograph;
	
	
	public MicrographParticle(int x, int y, Micrograph micrograph)
	{
		super(x, y);
		this.micrograph = micrograph;
	}
	
	
	
	public Micrograph getMicrograph() {
		return micrograph;
	}


	public void setMicrograph(Micrograph micrograph) {
		this.micrograph = micrograph;
	}


	
	


}
