package tiltpairpicker.model;

import java.util.ArrayList;
import java.util.List;

import trainingpicker.model.Micrograph;
import trainingpicker.model.Particle;

public class UntiltedMicrograph extends Micrograph {

	private TiltedMicrograph tiltedmicrograph;
	private List<UntiltedParticle> particles;

	public UntiltedMicrograph(String file, TiltedMicrograph tiltedmicrograph) {
		super(file);
		this.tiltedmicrograph = tiltedmicrograph;
		particles = new ArrayList<UntiltedParticle>();
	}

	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	public TiltedMicrograph getTiltedMicrograph() {
		return tiltedmicrograph;
	}

	public UntiltedParticle getParticle(int x, int y, int size) {
		for (UntiltedParticle p : particles)
			if (p.contains(x, y, size))
				return p;
		return null;

	}

	@Override
	public boolean hasData() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void reset() {
		// TODO Auto-generated method stub

	}

	public void removeParticle(UntiltedParticle p) {
		particles.remove(p);
		
	}

	public void addParticle(UntiltedParticle p) {
		particles.add(p);
		
	}
	
	public List<UntiltedParticle> getParticles()
	{
		return particles;
	}

}
