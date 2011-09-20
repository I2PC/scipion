package tiltpairpicker.model;

import java.util.ArrayList;
import java.util.List;

import trainingpicker.model.Micrograph;
import trainingpicker.model.Particle;
import xmipp.TiltPairAligner;

public class UntiltedMicrograph extends Micrograph {

	private TiltedMicrograph tiltedmicrograph;
	private List<UntiltedParticle> particles;
	private UntiltedParticle activeparticle;

	public UntiltedMicrograph(String file, TiltedMicrograph tiltedmicrograph) {
		super(file, getName(file, 1));
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
		tiltedmicrograph.removeParticle(p.getTiltedParticle());
		if(!particles.isEmpty())
			activeparticle = particles.get(particles.size() - 1);
		else
			activeparticle = null;
	}

	public void addParticle(UntiltedParticle p) {
		particles.add(p);
		activeparticle = p;
	}
	
	public List<UntiltedParticle> getParticles()
	{
		return particles;
	}
	
	public void setActiveParticle(UntiltedParticle p)
	{
		this.activeparticle = p;
	}

	public UntiltedParticle getActiveParticle()
	{
		return activeparticle;
	}

	public boolean hasActiveParticle()
	{
		return activeparticle != null;
	}

	public TiltedParticle getActiveTiltedParticle()
	{
		if(getActiveParticle() == null)
			return null;
		return activeparticle.getTiltedParticle();
			
	}
	
	

}
