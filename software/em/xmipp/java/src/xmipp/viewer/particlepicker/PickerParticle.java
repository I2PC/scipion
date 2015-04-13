package xmipp.viewer.particlepicker;

import xmipp.jni.Particle;

public class PickerParticle extends Particle
{

	private ParticleCanvas canvas;
	protected Micrograph micrograph;
        
        public PickerParticle(int x, int y, Micrograph m)
        {
            super(x, y);
            this.micrograph = m;
		// TODO Auto-generated constructor stub
        }
	
	public PickerParticle(int x, int y, Micrograph m, double cost)
	{
		super(x, y, cost);
		this.micrograph = m;
		// TODO Auto-generated constructor stub
	}
	
	public int getX0(int size)
	{
		int radius = size/2;
		return getX() - radius;
	}
	
	public int getY0(int size)
	{
		int radius = size/2;
		return getY() - radius;
	}
	
	
	public ParticleCanvas getParticleCanvas(ParticlePickerJFrame frame)
	{
		if(canvas == null)
			canvas = new ParticleCanvas(this, frame);
		return canvas; 
	}
	
	public void resetParticleCanvas()
	{
		canvas = null;
	}
	

	public Micrograph getMicrograph() {
		return micrograph;
	}


	public void setMicrograph(Micrograph micrograph) {
		this.micrograph = micrograph;
	}

}
